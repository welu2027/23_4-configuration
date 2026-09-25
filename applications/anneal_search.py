"""Simulated-annealing search for a combinatorial (23_4) with small sigma_2.

Uses the same edge switches as sample_random.py (about two per step) and accepts worse
moves with probability exp(-delta/T), with T decreasing linearly from 0.02. Seed 3 with
6,000 steps finds sigma_2 = 2.5633 < 2.5698 = sigma_2(B), the value quoted in Section 8.
Saves the incidence matrix of the best configuration to data/anneal_best.json.
"""
import json
import math
import os
import random
import sys

from common import load, switch_chain, sigma2, check_configuration


def anneal(seed, steps):
    rng = random.Random(seed)
    H = switch_chain(load('A'), rng, 20)
    cur = sigma2(H)
    best, bestH = cur, H.copy()
    for t in range(steps):
        T = 0.02 * (1 - t / steps) + 1e-4
        H2 = switch_chain(H, rng, 0.02)
        v = sigma2(H2)
        if v < cur or rng.random() < math.exp(-(v - cur) / T):
            H, cur = H2, v
            if cur < best:
                best, bestH = cur, H.copy()
    return best, bestH


if __name__ == '__main__':
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    steps = int(sys.argv[2]) if len(sys.argv) > 2 else 6000
    best, H = anneal(seed, steps)
    assert check_configuration(H)
    print(f'seed {seed}, {steps} steps: best sigma_2 = {best:.4f} (B has {sigma2(load("B")):.4f})')
    os.makedirs('data', exist_ok=True)
    with open('data/anneal_best.json', 'w') as f:
        json.dump({'seed': seed, 'steps': steps, 'sigma2': best, 'incidence': H.tolist()}, f)
