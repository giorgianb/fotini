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
        item = np.array(sum(p, ()))
        all_perms.append(item)


    return all_perms


