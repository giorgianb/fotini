import fotini as fi
from math import sqrt
import numpy as np
import tqdm
import matplotlib.pyplot as plt
from collections import defaultdict
from icecream import ic

bs = fi.BeamSplitter(1/sqrt(2), 1j/sqrt(2))
nodes = {0:bs, 1:bs}
edges = {0:{0:(1, 1, 1e-6*fi.Setup.c), 1:(1, 0, 2e-6*fi.Setup.c)}}
setup = fi.Setup(nodes=nodes, edges=edges)

φ_1 = fi.GaussianPhoton(ω_0=1e9, Δ=1, t_0=0, n=1)
probs = defaultdict(list)
ts = np.linspace(-10, 10, 250)

norms = []
for t in tqdm.tqdm(ts):
    φ_2 = fi.GaussianPhoton(ω_0=1e9, Δ=1, t_0=t, n=1)
    photon_inputs = [(φ_1, (0, 0)), (φ_2, (0, 1))]
    cp = setup.calculate(photon_inputs)

    for k, v in cp.items():
        probs[k].append(v)
    probs['norm'].append(sum(cp.values()))

for k, p in probs.items():
    plt.plot(ts, p, label=f'{k}')

plt.title('Mach-Zender Interferometer')
plt.ylim(0, 1.1)
plt.xlabel('$\\Delta t$')
plt.ylabel('$P(|x〉)$')
plt.legend()
plt.show()
