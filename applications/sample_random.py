"""Sample 20,000 random combinatorial (23_4) configurations and compare A and B.

Eight independent chains (seeds 1-8) of 2,500 samples each. Each chain starts from A,
runs 20 sweeps of burn-in, then records a sample every 3 sweeps. Output:
data/random_sample.json with sigma_2 and the minimum relay count of every sample.

Reproduces the numbers in Section 8: sigma_2 ranges from 2.621 to 3.028 in the sample,
the E-efficiency factor is at most 0.571, and the minimum relay count is 3, 4 or 5 in
21.0%, 76.6% and 2.4% of the sample (never 6).

Runtime: a few minutes on 8 cores.
"""
import json
import os
import random
from collections import Counter
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from common import load, switch_chain, sigma2, min_relays, e_efficiency, check_configuration

SEEDS = range(1, 9)
PER_SEED = 2500


def run(seed):
    rng = random.Random(seed)
    H = switch_chain(load('A'), rng, 20)
    out = []
    for _ in range(PER_SEED):
        H = switch_chain(H, rng, 3)
        assert check_configuration(H)
        out.append((sigma2(H), min_relays(H)))
    return out


if __name__ == '__main__':
    with ProcessPoolExecutor(8) as ex:
        res = [x for chunk in ex.map(run, SEEDS) for x in chunk]
    os.makedirs('data', exist_ok=True)
    with open('data/random_sample.json', 'w') as f:
        json.dump([{'sigma2': s, 'min_relays': r} for s, r in res], f)

    s = np.array([x[0] for x in res])
    rel = Counter(x[1] for x in res)
    N = len(res)
    print(f'samples: {N}')
    print(f'sigma_2 in sample: min {s.min():.4f}, max {s.max():.4f}')
    print(f'E-efficiency in sample: max {(16 - s.min() ** 2) / 16:.4f}')
    print('min relays: ' + ', '.join(f'{k}: {100 * v / N:.1f}%' for k, v in sorted(rel.items())))
    for name in 'AB':
        H = load(name)
        print(f'{name}: sigma_2 = {sigma2(H):.4f}, E-efficiency = {e_efficiency(H):.4f}, '
              f'min relays = {min_relays(H)}, sample configurations with smaller sigma_2: '
              f'{int((s < sigma2(H) - 1e-9).sum())}')
