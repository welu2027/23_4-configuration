from itertools import product, combinations
P = [(1, 0, 0)] + [(0, 1, s) for s in (1, -1)] + [(0, 3, 2 * s) for s in (1, -1)] + [(1, s, 0) for s in (1, -1)]\
    + [(a, e * b, d * c) for (a, b, c) in [(3, 2, 2), (2, 3, 2), (6, 1, 2), (6, 5, 2)]
       for e, d in product((1, -1), repeat=2)]
L = [(1, 0, 0)] + [(0, 1, s) for s in (1, -1)] + [(0, 2, 3 * s) for s in (1, -1)] + [(1, 0, 3 * s) for s in (1, -1)]\
    + [(a, e * b, d * c) for (a, b, c) in [(1, 2, 2), (2, 6, 9), (2, 2, 1), (2, 2, 5)]
       for e, d in product((1, -1), repeat=2)]
assert len(P) == 23 and len(L) == 23
inc = [[sum(x * y for x, y in zip(p, l)) == 0 for l in L] for p in P]
assert all(sum(row) == 4 for row in inc)
assert all(sum(inc[i][j] for i in range(23)) == 4 for j in range(23))
assert all(sum(a and b for a, b in zip(inc[i], inc[j])) <= 1 for i, j in combinations(range(23), 2))
cross = lambda u, v: (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])
assert all(cross(P[i], P[j]) != (0, 0, 0) for i, j in combinations(range(23), 2))
assert all(cross(L[i], L[j]) != (0, 0, 0) for i, j in combinations(range(23), 2))
print("geometric (23_4) configuration verified")
