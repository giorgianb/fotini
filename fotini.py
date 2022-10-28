from sympy import MatrixSymbol
from sympy import symbols, simplify, expand
import sympy
import numpy as np
from collections import defaultdict, namedtuple
from scipy.integrate import nquad
from itertools import permutations, product
import math
import time
from icecream import ic
import ctypes
import os
from scipy import LowLevelCallable
import utils

class Data(ctypes.Structure):
    _fields_ = [
            ('n_terms', ctypes.c_int),
            ('norms', ctypes.POINTER(ctypes.c_double)),
            ('coeffs_real', ctypes.POINTER(ctypes.c_double)),
            ('coeffs_imag', ctypes.POINTER(ctypes.c_double)),
            ('n_permutation', ctypes.c_int),
            ('permutations', ctypes.POINTER(ctypes.c_int)),
            ('t_0', ctypes.POINTER(ctypes.c_double)),
            ('ω_0', ctypes.POINTER(ctypes.c_double)),
            ('Δ', ctypes.POINTER(ctypes.c_double)),
            ('printed', ctypes.c_int),
    ]


GaussianPhoton = namedtuple('GaussianPhoton', 't_0 ω_0 Δ')
BeamSplitter = namedtuple('BeamSplitter', 'R T')

class Setup:
    def __init__(self, nodes, edges):
        self.nodes = nodes
        self.edges = edges

        self.input_operators = {}
        self.output_operators = {}
        self.input_indices = {}
        self.output_indices = {}
        for node in nodes.keys():
            # Each beamsplitter has two inputs and two outputs
            for i in range(2):
                input_name = f'a_i_{node}_{i}'
                output_name = f'a_o_{node}_{i}'
                input_symbol = symbols(input_name, commutative=False)
                output_symbol = symbols(output_name, commutative=False)
                self.input_indices[(node, i)] = input_symbol
                self.input_operators[input_symbol] = (node, i)
                self.output_indices[(node, i)] = output_symbol
                self.output_operators[output_symbol] = (node, i)

        self.expanded_operators = {}
        self.outputs = set()
        dependent = set()
        def resolve(operator, visited=set()):
            #print(f'resolve(operator={operator}, visited={visited})')
            if operator in visited:
                raise ValueError("Edges form a cycle")

            if operator in self.expanded_operators:
                return self.expanded_operators[operator]
            elif operator in self.input_operators:
                if len(visited) != 0:
                    dependent.add(self.input_operators[operator])

                node, input_index = self.input_operators[operator]
                bs = self.nodes[node]
                a_o0 = resolve(self.output_indices[(node, 0)], visited | {operator})
                a_o1 = resolve(self.output_indices[(node, 1)], visited | {operator})
                if input_index == 0:
                    a_cur = bs.T * a_o0 + bs.R * a_o1
                else:
                    a_cur = bs.T * a_o1 + bs.R * a_o0
                self.expanded_operators[operator] = a_cur
                return a_cur
            elif operator in self.output_operators:
                node, output_index = self.output_operators[operator]
                if node not in self.edges:
                    self.expanded_operators[operator] = operator
                    self.outputs.add((node, output_index))
                    return operator

                connections = self.edges[node]
                if output_index not in connections:
                    self.expanded_operators[operator] = operator
                    self.outputs.add((node, output_index))

                connection_node, input_index = connections[output_index]
                a_in = resolve(self.input_indices[(connection_node, input_index)], visited | {operator})
                self.expanded_operators[operator] = a_in
                return a_in

        for operator in self.input_operators.keys():
            resolve(operator, visited=set())
        self.inputs = set(self.input_indices.keys()) - dependent

    def calculate(self, photon_inputs):
        expr = 1
        ω_0 = []
        Δ = []
        t_0 = []
        for i, (photon, input_index) in enumerate(photon_inputs):
            if input_index not in self.inputs:
                raise ValueError("Attempting to create a photon at a non-input node")
            ω_0.append(photon.ω_0)
            Δ.append(photon.Δ)
            t_0.append(photon.t_0)

            expr *= self.expanded_operators[self.input_indices[input_index]]
            t = symbols(f't_{i}', commutative=False)
            expr *= t

        expr = expand(expr)
        ω_0 = np.array(ω_0, dtype=ctypes.c_double)
        Δ = np.array(Δ, dtype=ctypes.c_double)
        t_0 = np.array(t_0, dtype=ctypes.c_double)

        count_mapping = defaultdict(list)
        output_mapping = {v:i for i, v in enumerate(self.outputs)}
        for term in expr.args:
            coefficient = term.args[0]
            order = []
            count = [0]*len(output_mapping)
            for subterm in term.args[1:]:
                if subterm in self.output_operators:
                    output_index = self.output_operators[subterm]
                    count[output_mapping[output_index]] += 1
                    # Maybe check if this is right? Because we do define the output order
                    # somewhat arbitrarily
                    order.append(output_mapping[output_index])
            count_mapping[tuple(count)].append((coefficient, tuple(order)))


        count_mapping = dict(count_mapping)
        def generate_term(terms):
            coefficients = []
            norms = []
            permutations = []
            for coefficient, orig_order in terms:
                coefficients.append(coefficient)
                norms.append(utils.get_normalization_coefficient(orig_order))
                permutations.append(utils.generate_permutation(orig_order))

            return coefficients, norms, permutations

        count_probabilities = {}
        data = Data()
        user_data = ctypes.cast(ctypes.pointer(data), ctypes.c_void_p)
        data.t_0 = (ctypes.c_double*len(self.inputs)).from_buffer(t_0);
        data.ω_0 = (ctypes.c_double*len(self.inputs)).from_buffer(ω_0);
        data.Δ = (ctypes.c_double*len(self.inputs)).from_buffer(Δ);

        lib = ctypes.CDLL(os.path.abspath('integrands.so'))
        lib.β.restype = ctypes.c_double
        lib.β.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.c_void_p)
        func = LowLevelCallable(lib.β, user_data)
        for count, terms in count_mapping.items():
            coefficients, norms, permutations = generate_term(terms)
            print(f'[Integrating for {count}]')
            n_permutation = len(permutations[0])
            coeffs_real = np.array(coefficients).real.astype(ctypes.c_double)
            coeffs_imag = np.array(coefficients).imag.astype(ctypes.c_double)
            norms = np.array(norms, dtype=ctypes.c_double)
            permutations = np.stack(permutations).astype(ctypes.c_int).reshape(-1, 4)

            data.n_terms = ctypes.c_int(len(terms));
            data.norms = (ctypes.c_double * len(terms)).from_buffer(norms)
            data.coeffs_real = (ctypes.c_double * len(terms)).from_buffer(coeffs_real)
            data.coeffs_imag = (ctypes.c_double * len(terms)).from_buffer(coeffs_imag)
            data.permutations = (ctypes.c_int*len(permutations)).from_buffer(permutations)
            data.n_permutation = n_permutation
            data.printed = ctypes.c_int(0)

            t1 = time.time()
            count_probabilities[count] = nquad(
                    func, 
                    [(-np.inf, np.inf)]*4, 
            )[0]
            t2 = time.time()
            print(f'Probability = {count_probabilities[count]}\nTime = {t2 - t1:3f}s')

        return count_probabilities
