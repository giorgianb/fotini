# Some boilerplate so the imports work properly
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), "../"))

import fotini as fi
from math import sqrt
import numpy as np
import tqdm
import matplotlib.pyplot as plt
from collections import defaultdict
from icecream import ic

N_SAMPLES = 100
bs = fi.BeamSplitter(1/sqrt(2), 1j/sqrt(2))
nodes = {0:bs, 1:bs}

probs = defaultdict(list)
distances = np.linspace(0, 6e-7, N_SAMPLES)
for distance in tqdm.tqdm(distances):
    # Step 1: Specify the setup. Here we vary the length of the second arm
    edges = {0:{0:(1, 1, 1e-6*fi.Setup.c), 1:(1, 0, distance*fi.Setup.c)}}
    setup = fi.Setup(nodes=nodes, edges=edges)

    # Step 2: place the photon at the input
    φ = fi.GaussianPhoton(ω_0=1e9, Δ=1, t_0=0, n=1)
    photon_inputs = [(φ, (0, 0))]

    # Step 3: calculate the probabilities
    cp = setup.calculate(photon_inputs)

    for k, v in cp.items():
        probs[k].append(v)

plt.figure(figsize=(10, 5))
plt.plot(distances, probs[(0, 1)], 'r', label=f'(0, 1)')
plt.plot(distances, probs[(1, 0)], 'k', label=f'(1, 0)')

plt.title('Mach-Zender Interferometer')
plt.ylim(0, 1.1)
plt.xlabel('$d$')
plt.ylabel('$P(|x〉)$')
plt.legend()
plt.savefig('figs/mz-0.png')
