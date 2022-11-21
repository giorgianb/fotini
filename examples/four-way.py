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
nodes = {0:bs, 1:bs, 2:bs}
edges = {
        # Format: (output_port: (input_node, input_port, distance))
        # Specify the edges for node 0
        0: {0: (2, 0, 1), 1: (1, 1, 1)},
}

setup = fi.Setup(nodes=nodes, edges=edges)

# Step 2: Specify Inputs
φ_1 = fi.GaussianPhoton(ω_0=1e3, Δ=1, t_0=1, n=1)
φ_2 = fi.GaussianPhoton(ω_0=1e3, Δ=1, t_0=0, n=1)
φ_4 = fi.GaussianPhoton(ω_0=1e3, Δ=1, t_0=1, n=1)
probs = defaultdict(list)
ts = np.linspace(-5, 5, 250)


norms = []
for t in tqdm.tqdm(ts):
    φ_3 = fi.GaussianPhoton(ω_0=1e3, Δ=1, t_0=t, n=1)
    # Specify the photons that are going to each input port
    photon_inputs = [(φ_1, (1, 0)), (φ_2, (0, 0)), (φ_3, (0, 1)), (φ_4, (2, 1)),]

    # Step 3: Compute output staistics for a specific set of inputs
    cp = setup.calculate(photon_inputs)

    for k, v in cp.items():
        probs[k].append(v)
    probs['norm'].append(sum(cp.values()))


plt.figure(figsize=(10, 5))
for k, p in probs.items():
    plt.plot(ts, p, label=f'{k}')

plt.title('Four-Way Output Statistics')
plt.ylim(0, 0.1)
plt.xlabel('$\\Delta t$')
plt.ylabel('$P(|x〉)$')
plt.savefig('figs/four-way.png')
