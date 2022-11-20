from collections import defaultdict
import math
import numpy as np
from itertools import permutations, product

def get_normalization_coefficient(operators):
    counts = defaultdict(int)
    for item in operators:
        counts[item] += 1

    n = 1
    for item in counts.values():
        n *= math.sqrt(math.factorial(item))


    return n

def generate_permutation(order):
    positions = defaultdict(list)
    for i, t in enumerate(order):
        positions[t].append(i)

    perms = []
    for k, v in sorted(positions.items()):
        perms.append(permutations(v))

    all_perms = []
    for p in product(*perms):
        order = np.array(sum(p, ()))
        perm = np.zeros(len(order), dtype=int)
        perm[order] = np.arange(len(order))
        all_perms.append(perm)


    return all_perms

def integrate(t_0, ω_0, Δ, coefficients, norms, permutations):
    matrix = np.stack([-Δ**2, 2*t_0*Δ**2 - 1j*ω_0, -t_0**2*Δ**2], axis=-1)
    coeffs = np.prod((2*Δ**2/np.pi)**(1/4))*coefficients*norms/permutations.shape[1]
    total = 0
    for c_0, perms_0 in zip(coeffs, permutations):
        for c_1, perms_1 in zip(coeffs, permutations):
            for p_0 in perms_0:
                for p_1 in perms_1:
                    total += c_0*c_1.conjugate()*integrate_wavelet(matrix[p_0] + matrix[p_1].conjugate())

    return total.real

def integrate_wavelet(matrix):
    total = np.sum(matrix[:, 2] - matrix[:, 1]**2/(4*matrix[:, 0]))
    prod = np.prod(np.sqrt(np.pi)/np.sqrt(-matrix[:, 0]))

    return prod*np.exp(total)
