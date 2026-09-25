# Connections to other fields (Section 8)

Code and data for Section 8 of the paper. It compares the two geometric $(23_4)$
configurations with random combinatorial $(23_4)$ configurations, as expander graphs,
as incomplete block designs and as key predistribution schemes.

Run everything from this folder. The configurations are read from
`../integer_23_4_certificate.json` (A) and `../geometric_23_4_certificate.json` (B).

Requirements: Python 3, `numpy`, `networkx`, `sympy`, `matplotlib`.

## Scripts

| Script | What it does | Runtime |
|---|---|---|
| `exact_B.py` | Proves B's two properties exactly: the factored characteristic polynomial (σ₂² is the largest root of y³ − 11y² + 31y − 13), and the 18 orbit representatives of non-collinear pairs under the order-8 automorphism group, each with 6, 7 or 8 relays. | seconds |
| `sample_random.py` | Samples 20,000 random combinatorial $(23_4)$ configurations by degree- and girth-preserving edge switches of the Levi graph (8 seeded chains of 2,500). Writes `data/random_sample.json`. | a few minutes on 8 cores |
| `anneal_search.py` | Simulated-annealing search for small σ₂. Seed 3 with 6,000 steps finds σ₂ = 2.5633. Writes `data/anneal_best.json`. | seconds |
| `make_figure.py` | Draws the two-panel comparison figure from the data files. | seconds |
| `common.py` | Shared helpers: loading, edge switches, σ₂, E-efficiency, relay counts. | |

`data/` contains the output of `sample_random.py` and `anneal_search.py`, so the figure
and the statistics can be checked without rerunning the sampling.

## Results

| | σ₂ | E-efficiency | fewest relays |
|---|---|---|---|
| B (self-polar) | 2.5698 | 0.587 | 6 |
| A (integral) | 2.8853 | 0.480 | 3 |
| 20,000 random configurations | 2.621 to 3.028 | at most 0.571 | 3, 4 or 5 in 21.0%, 76.6%, 2.4% |
| best found by annealing | 2.5633 | 0.589 | 5 |

Every random configuration has larger σ₂ than B, and none has six relays. The random
configurations need not be geometric.

## Definitions

- **σ₂** is the second largest singular value of the 23 × 23 incidence matrix $N$ (lines by points). Smaller means the Levi graph is a better expander.
- **E-efficiency factor.** With the points as 23 treatments and the lines as 23 blocks of size four, $N^TN = 4I + M$, where $M$ is the collinearity adjacency matrix, so the factor is $(16-\sigma_2^2)/16$.
- **Relays.** With the points as nodes and the lines as keys, two nodes share a key exactly when they are collinear. A relay for two non-collinear nodes is a node collinear with both.
