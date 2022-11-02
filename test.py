import fotini as fi
from math import sqrt
import numpy as np
import tqdm
import matplotlib.pyplot as plt
from collections import defaultdict
from icecream import ic

bs = fi.BeamSplitter(1/sqrt(2), 1j/sqrt(2))
setup = fi.Setup(nodes={0:bs}, edges={})

φ_1 = fi.GaussianPhoton(ω_0=1e3, Δ=1, t_0=0, n=1)
probs = defaultdict(list)
ts = np.linspace(-10, 10, 250)



norms = []
for t in tqdm.tqdm(ts):
    φ_2 = fi.GaussianPhoton(ω_0=1e3, Δ=1, t_0=t, n=1)
    photon_inputs = [(φ_1, (0, 0)), (φ_2, (0, 1))]
    cp = setup.calculate(photon_inputs)

    for k, v in cp.items():
        probs[k].append(v)
    probs['norm'].append(sum(cp.values()))

for k, p in probs.items():
    plt.plot(ts, p, label=f'{k}')

plt.title('Hong-Ou Mandel Effect')
plt.ylim(0, 1.1)
plt.xlabel('$\\Delta t$')
plt.ylabel('$P(|x〉)$')
plt.legend()
plt.show()

