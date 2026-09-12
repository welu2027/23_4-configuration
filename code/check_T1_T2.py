import json, os, itertools
from fractions import Fraction as F

def t_mul(a, b):
    r = [F(0)] * 7
    for i in range(4):
        if a[i]:
            for j in range(4): r[i + j] += a[i] * b[j]
    for k in (6, 5, 4):
        if r[k]:
            c = r[k]; r[k] = F(0); r[k - 2] += c; r[k - 4] += c
    return tuple(r[:4])
def t_add(a, b): return tuple(x + y for x, y in zip(a, b))
def t_sub(a, b): return tuple(x - y for x, y in zip(a, b))
def t_zero(a): return all(x == 0 for x in a)
T = lambda *c: tuple(F(x) for x in c) + (F(0),) * (4 - len(c))

def q_mul(a, b):
    a0, a1, a2, a3 = a; b0, b1, b2, b3 = b
    return (a0*b0 + 2*a1*b1 + 5*a2*b2 + 10*a3*b3,
            a0*b1 + a1*b0 + 5*a2*b3 + 5*a3*b2,
            a0*b2 + a2*b0 + 2*a1*b3 + 2*a3*b1,
            a0*b3 + a3*b0 + a1*b2 + a2*b1)
def q_add(a, b): return tuple(x + y for x, y in zip(a, b))
def q_sub(a, b): return tuple(x - y for x, y in zip(a, b))
def q_zero(a): return all(x == 0 for x in a)
Q = lambda a=0, b=0, c=0, d=0: (F(a), F(b), F(c), F(d))

def run(name, P0, L0, mul, add, sub, zero, one, minus_one):
    V = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
    sgn = lambda s, x: x if s == 1 else tuple(-c for c in x)
    def cross(u, v):
        return (sub(mul(u[1], v[2]), mul(u[2], v[1])),
                sub(mul(u[2], v[0]), mul(u[0], v[2])),
                sub(mul(u[0], v[1]), mul(u[1], v[0])))
    same = lambda p, q: all(zero(c) for c in cross(p, q))
    def orbit(p):
        o = []
        for v in V:
            q = tuple(sgn(s, c) for s, c in zip(v, p))
            if not any(same(q, x) for x in o): o.append(q)
        return o
    P, L, sp, sl = [], [], [], []
    for p in P0: o = orbit(p); sp.append(len(o)); P += o
    for l in L0: o = orbit(l); sl.append(len(o)); L += o
    dot = lambda p, l: add(add(mul(p[0], l[0]), mul(p[1], l[1])), mul(p[2], l[2]))
    inc = [[zero(dot(p, l)) for l in L] for p in P]
    print("%s: point orbit sizes %s -> %d points; line orbit sizes %s -> %d lines"
          % (name, sp, len(P), sl, len(L)))
    print("   22 distinct points: %s   22 distinct lines: %s"
          % (all(not same(P[i], P[j]) for i, j in itertools.combinations(range(len(P)), 2)),
             all(not same(L[i], L[j]) for i, j in itertools.combinations(range(len(L)), 2))))
    print("   every point on exactly 4 lines: %s" % all(sum(r) == 4 for r in inc))
    print("   every line has exactly 4 points: %s"
          % all(sum(inc[i][j] for i in range(len(P))) == 4 for j in range(len(L))))
    print("   no two points on two common lines: %s"
          % all(sum(a and b for a, b in zip(inc[i], inc[j])) <= 1
                for i, j in itertools.combinations(range(len(P)), 2)))
    return inc

P1 = [(T(0), T(1), T(2, 0, -1)), (T(1), T(0), T(-1, 0, 1)), (T(1), T(0, 0, -1), T(0)),
      (T(1), T(1), T(1)), (T(1), T(0, 1, 1), T(0, 1, 1)), (T(1), T(0, 0, 0, 1), T(1)),
      (T(1), T(1), T(1, 1, -1, -1))]
L1 = [(T(0), T(1), T(-1)), (T(1), T(0), T(-1)), (T(1), T(-1), T(0)),
      (T(1), T(-1, 0, 1), T(0, 0, -1)), (T(1), T(1, 1, 0, -1), T(-1, 0, -1, 1)),
      (T(1), T(0, 1), T(0, 0, -1)), (T(1), T(-1, 0, 1), T(1, -1))]
inc1 = run("T_1 over Q(sqrt(phi))", P1, L1, t_mul, t_add, t_sub, t_zero, T(1), T(-1))

h = F(1, 2)
P2 = [(Q(), Q(1), Q(-h, 0, h)), (Q(1), Q(), Q(F(3,2), 0, -h)), (Q(1), Q(h, 0, -h), Q()),
      (Q(1), Q(1), Q(1)), (Q(1), Q(2, F(3,2), -1, -h), Q(2, F(3,2), -1, -h)),
      (Q(1), Q(1, 1), Q(1)), (Q(1), Q(1), Q(1, h, 0, -h))]
L2 = [(Q(), Q(1), Q(-1)), (Q(1), Q(), Q(-1)), (Q(1), Q(-1), Q()),
      (Q(1), Q(h, 0, h), Q(-F(3,2), 0, -h)),
      (Q(1), Q(F(3,2), -h, -h, h), Q(h, -1, -h)),
      (Q(1), Q(h, -h, h, -h), Q(-F(3,2), 0, -h)),
      (Q(1), Q(h, 0, h), Q(-F(3,2), -h, -h, -h))]
inc2 = run("T_2 over Q(sqrt2,sqrt5)", P2, L2, q_mul, q_add, q_sub, q_zero, Q(1), Q(-1))

try:
    import networkx as nx
    from networkx.algorithms.isomorphism import GraphMatcher
    def levi(m):
        G = nx.Graph(); k = len(m)
        for i in range(k): G.add_node(("p", i), side=0); G.add_node(("l", i), side=1)
        for i in range(k):
            for j in range(k):
                if m[i][j]: G.add_edge(("p", i), ("l", j))
        return G
    GC = None
    for fn in ("cuntz22_inc.json", "../certificates/cuntz22_inc.json"):
        if os.path.exists(fn): GC = levi(json.load(open(fn))); break
    sides = lambda a, b: a["side"] == b["side"]
    for nm, m in (("T_1", inc1), ("T_2", inc2)):
        G = levi(m)
        print("%s: Levi girth %d, Aut kept apart %d, with dualities %d"
              % (nm, nx.girth(G),
                 sum(1 for _ in GraphMatcher(G, G, node_match=sides).isomorphisms_iter()),
                 sum(1 for _ in GraphMatcher(G, G).isomorphisms_iter())))
        if GC is not None:
            print("   isomorphic to Cuntz's (22_4): %s" % nx.is_isomorphic(G, GC))
    if GC is None: print("(cuntz22_inc.json not found; run check_cuntz.py for the comparison)")
except ImportError:
    print("(networkx not installed; skipping Levi automorphism computation)")
