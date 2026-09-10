import json, subprocess, re, itertools, math, mpmath, pynauty
from fractions import Fraction as F
from classify_Vn import Layout, build, variables, poly_str, singular_script, sym_coords
mpmath.mp.dps = 40
class K:
    @staticmethod
    def mul(x, y): return (x[0] * y[0] + 17 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])
    @staticmethod
    def add(x, y): return (x[0] + y[0], x[1] + y[1])
    @staticmethod
    def sub(x, y): return (x[0] - y[0], x[1] - y[1])
    @staticmethod
    def inv(x): n = x[0] ** 2 - 17 * x[1] ** 2; return (x[0] / n, -x[1] / n)
    @staticmethod
    def zero(x): return x == (0, 0)
ONE, ZERO = (F(1), F(0)), (F(0), F(0))
def cross(u, v):
    return (K.sub(K.mul(u[1], v[2]), K.mul(u[2], v[1])), K.sub(K.mul(u[2], v[0]), K.mul(u[0], v[2])), K.sub(K.mul(u[0], v[1]), K.mul(u[1], v[0])))
def matvec(M, p): return tuple(K.add(K.add(K.mul(M[i][0], p[0]), K.mul(M[i][1], p[1])), K.mul(M[i][2], p[2])) for i in range(3))
def solve3(A, b):
    M = [list(A[i]) + [b[i]] for i in range(3)]
    for c in range(3):
        piv = next(r for r in range(c, 3) if not K.zero(M[r][c])); M[c], M[piv] = M[piv], M[c]
        inv = K.inv(M[c][c]); M[c] = [K.mul(x, inv) for x in M[c]]
        for r in range(3):
            if r != c and not K.zero(M[r][c]):
                f = M[r][c]; M[r] = [K.sub(x, K.mul(f, y)) for x, y in zip(M[r], M[c])]
    return [M[r][3] for r in range(3)]
def projectivity(src, dst):
    def frame(pts):
        A = [[pts[j][i] for j in range(3)] for i in range(3)]
        a = solve3(A, pts[3]); return [[K.mul(a[j], pts[j][i]) for j in range(3)] for i in range(3)]
    S, D = frame(src), frame(dst)
    cols = [solve3(S, [ONE if i == k else ZERO for i in range(3)]) for k in range(3)]
    Sinv = [[cols[k][i] for k in range(3)] for i in range(3)]
    return [[K.add(K.add(K.mul(D[i][0], Sinv[0][j]), K.mul(D[i][1], Sinv[1][j])), K.mul(D[i][2], Sinv[2][j])) for j in range(3)] for i in range(3)]
def z(a, b): return (F(a), F(b))
B = [(z(1,0), z(0,0), z(0,0)), (z(0,0), z(4,0), z(4,0)), (z(0,0), z(4,0), z(-4,0)), (z(0,0), z(4,0), z(-1,1)), (z(0,0), z(4,0), z(1,-1)), (z(8,0), z(4,0), z(0,0)), (z(-8,0), z(4,0), z(0,0))]
for (X, Z) in ((z(4,0), z(3,1)), (z(2,2), z(4,0)), (z(4,0), z(5,-1)), (z(-2,2), z(-1,1))):
    for e, d in itertools.product((1, -1), repeat=2): B.append(((e * X[0], e * X[1]), z(4,0), (d * Z[0], d * Z[1])))
polarB = lambda p: ((2 * p[0][0], 2 * p[0][1]), (-4 * p[1][0], -4 * p[1][1]), tuple(-x for x in K.mul(z(1,1), p[2])))
LB = [polarB(p) for p in B]
dot = lambda p, l: K.add(K.add(K.mul(p[0], l[0]), K.mul(p[1], l[1])), K.mul(p[2], l[2]))
incB = [[int(K.zero(dot(p, l))) for l in LB] for p in B]
assert all(sum(r) == 4 for r in incB)
data = json.load(open('kn/23_B_4_4_3_3.json')); cs = data['case']; lay = Layout(23, cs[0], cs[1], cs[2], cs[3], cs[4]); vs = variables(lay)
key = sorted(data['structures'])[119]; eqs, _ = build(lay, data['structures'][key])
script = 'LIB "solve.lib";\n' + singular_script(lay, eqs, 'X') + 'quit;\n'; open('/tmp/e119.sing', 'w').write(script)
out = subprocess.run(['Singular', '-q', '/tmp/e119.sing'], capture_output=True, text=True, timeout=60).stdout
blk = re.search(r'SOLS X\n(.*?)ENDSOLS', out, flags=re.S).group(1)
vals = [l.strip() for l in blk.splitlines() if l.strip() and not re.match(r'^\[\d+\]:$', l.strip())]
sols = [[mpmath.mpf(vals[k + j].replace(' ', '')) for j in range(len(vs))] for k in range(0, len(vals) - len(vs) + 1, len(vs))]
def ident(x):
    r = mpmath.pslq([x, 1, mpmath.sqrt(17)], maxcoeff=10 ** 6, maxsteps=10 ** 5); a, b = F(-r[1], r[0]), F(-r[2], r[0])
    assert abs(mpmath.mpf(a.numerator) / a.denominator + mpmath.mpf(b.numerator) / b.denominator * mpmath.sqrt(17) - x) < mpmath.mpf(10) ** -30
    return (a, b)
pts_s, lns_s = sym_coords(lay)
def levi_label(inc):
    g = pynauty.Graph(46, directed=False, vertex_coloring=[set(range(23)), set(range(23, 46))]); g.set_adjacency_dict({p: [23 + l for l in range(23) if inc[p][l]] for p in range(23)})
    gens, o1, o2, _, _ = pynauty.autgrp(g); return pynauty.canon_label(g), gens
labB, gensB = levi_label(incB)
def group(gens):
    els = {tuple(range(46))}; fr = [tuple(range(46))]
    while fr:
        new = []
        for g in fr:
            for h in gens:
                x = tuple(h[g[i]] for i in range(46))
                if x not in els: els.add(x); new.append(x)
        fr = new
    return els
GB = group(gensB)
report = []
for si, s in enumerate(sols):
    ex = {vs[j]: ident(s[j]) for j in range(len(vs))}; ex['u1'] = ONE; ex['v1'] = ONE
    def ev(c):
        sg, var = c
        if sg == 0: return ZERO
        x = ex[var] if var else ONE; return (sg * x[0], sg * x[1])
    P = [tuple(ev(c) for c in pts_s[i]) for i in range(23)]; L = [tuple(ev(c) for c in lns_s[i]) for i in range(23)]
    inc = [[int(K.zero(dot(p, l))) for l in L] for p in P]
    assert all(sum(r) == 4 for r in inc) and all(sum(inc[i][j] for i in range(23)) == 4 for j in range(23))
    lab, _ = levi_label(inc); invB = {v: i for i, v in enumerate(labB)}
    base = [labB.index(lab[i]) if False else None for i in range(23)]
    pos_s = {v: k for k, v in enumerate(lab)}
    base = [labB[pos_s[i]] for i in range(23)]
    found = None
    for g in GB:
        perm = [g[base[i]] for i in range(23)]
        if any(perm[i] >= 23 for i in range(23)): continue
        quad = None
        for q in itertools.combinations(range(23), 4):
            A = [P[i] for i in q]
            if all(not K.zero(dot(cross(A[a], A[b]), A[c])) for a, b, c in itertools.combinations(range(4), 3)): quad = q; break
        try: M = projectivity([P[i] for i in quad], [B[perm[i]] for i in quad])
        except StopIteration: continue
        if all(all(K.zero(x) for x in cross(matvec(M, P[i]), B[perm[i]])) for i in range(23)): found = (M, perm); break
    ok = found is not None
    print("solution %d (t1 = %s + %s sqrt17): exact projectivity to (B) found: %s" % (si, ex['t1'][0], ex['t1'][1], ok))
    if ok:
        M, perm = found; scale = K.inv(next(x for row in M for x in row if not K.zero(x)))
        Mn = [[K.mul(x, scale) for x in row] for row in M]
        print("   M =", [["%s+%s*s17" % (x[0], x[1]) for x in row] for row in Mn])
        report.append({'solution': si, 'coordinates': {k: [str(v[0]), str(v[1])] for k, v in ex.items()}, 'matrix': [[[str(x[0]), str(x[1])] for x in row] for row in Mn], 'point_map': perm})
json.dump({'note': 'M maps the 23 points of the structure-119 solution (in the normalization u1=v1=1, coordinates a+b*sqrt17) onto the 23 points of configuration (B) in its Z[sqrt17] model; verified exactly on all 23 points', 'solutions': report}, open('/home/claude/exact_equivalence_119.json', 'w'))
