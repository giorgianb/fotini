# Some boilerplate so the imports work properly
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), "../"))

import fotini as fi
from math import sqrt
import numpy as np
import tqdm
import matplotlib.pyplot as plt
from collections import defaultdict

# Step 1: Specify Set-Up
bs = fi.BeamSplitter(1/sqrt(2), 1j/sqrt(2))
setup = fi.Setup(nodes={0:bs}, edges={})

# Step 2: Specify Inputs
φ_1 = fi.GaussianPhoton(ω_0=1e3, Δ=1, t_0=0, n=1)
probs = defaultdict(list)
ts = np.linspace(-5, 5, 250)


norms = []
for t in tqdm.tqdm(ts):
    φ_2 = fi.GaussianPhoton(ω_0=1e3, Δ=1, t_0=t, n=1)
    # Specify the photons that are going to each input port
    photon_inputs = [(φ_1, (0, 0)), (φ_2, (0, 1))]

    # Step 3: Compute output staistics for a specific set of inputs
    cp = setup.calculate(photon_inputs)

    for k, v in cp.items():
        probs[k].append(v)
    probs['norm'].append(sum(cp.values()))

# Plot the computed output statistics
plt.figure(figsize=(10, 5))
plt.plot(ts, probs[(1, 1)], 'r', label='(1, 1)')
plt.plot(ts, probs[(0, 2)], 'b', linewidth=3, label='(0, 2)')
plt.plot(ts, probs[(2, 0)], 'r--', linewidth=2, label='(2, 0)')

plt.title('Hong-Ou Mandel Effect')
plt.ylim(0, 0.8)
plt.xlabel('$\\Delta t$')
plt.ylabel('$P(|x〉)$')
plt.legend()
plt.savefig('figs/hom-0.png')
