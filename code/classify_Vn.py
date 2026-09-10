import itertools, json, sys, time, re, subprocess, math
import pynauty
V = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
E, S1, S2, S3 = V
VIDX = {v: i for i, v in enumerate(V)}
def vmul(a, b): return tuple(x * y for x, y in zip(a, b))
SUB = {1: (E, S1), 2: (E, S2), 3: (E, S3)}
SIG = {1: S1, 2: S2, 3: S3}
OTHER = {1: S2, 2: S1, 3: S1}
class Layout:
    def __init__(self, n, kind, c, cp, bpairs, bplines):
        self.n, self.kind, self.c, self.cp = n, kind, c, cp
        self.centers = {'A': [], 'B': [1], 'C': [1, 2, 3]}[kind]
        self.tpairs = {'A': [], 'B': [1, 1], 'C': [1, 2, 3]}[kind]
        self.bpairs, self.bplines = list(bpairs), list(bplines)
        self.nc, self.nt, self.nb, self.nbl = len(self.centers), len(self.tpairs), len(self.bpairs), len(self.bplines)
        self.pstart = self.nc + 2 * (self.nt + self.nb); self.lstart = self.nc + 2 * (self.nt + self.nbl)
        assert self.pstart + 4 * c == n and self.lstart + 4 * cp == n, (n, kind, c, cp, bpairs, bplines)
        self.pair_axes = self.tpairs + self.bpairs
        self.line_centers = self.tpairs + self.bplines
    def O(self, i): return self.centers.index(i)
    def pair(self, k, pr): return self.nc + 2 * k + pr
    def P(self, j, a): return self.pstart + 4 * (j - 1) + VIDX[a]
    def axis(self, i): return self.centers.index(i)
    def pline(self, k, pr): return self.nc + 2 * k + pr
    def G(self, a, al): return self.lstart + 4 * (a - 1) + VIDX[al]
    def act_pt(self, s, v):
        if v < self.nc: return v
        if v < self.pstart:
            k, pr = divmod(v - self.nc, 2); return v if s in SUB[self.pair_axes[k]] else self.pair(k, 1 - pr)
        j = (v - self.pstart) // 4 + 1; al = V[(v - self.pstart) % 4]; return self.P(j, vmul(al, s))
    def act_ln(self, s, v):
        if v < self.nc: return v
        if v < self.lstart:
            k, pr = divmod(v - self.nc, 2); return v if s in SUB[self.line_centers[k]] else self.pline(k, 1 - pr)
        a = (v - self.lstart) // 4 + 1; al = V[(v - self.lstart) % 4]; return self.G(a, vmul(al, s))
def coset_reps(H):
    reps, seen = [], set()
    for v in V:
        if v in seen: continue
        reps.append(v); seen |= {vmul(v, h) for h in H}
    return reps
def attachments(n_targets, kinds, sizes, nt=0):
    states = {None: ()}
    for k in range(len(kinds)):
        nxt = {}
        for att in states.values():
            used = {(kinds[i], t) for i in range(len(att)) for t in att[i]}
            for pr in itertools.combinations(range(1, n_targets + 1), sizes[k]):
                if any((kinds[k], t) in used for t in pr): continue
                new = att + (pr,)
                N = n_targets + len(new); adj = {i: [] for i in range(N)}
                for i, p in enumerate(new):
                    for t in p: adj[n_targets + i].append(t - 1); adj[t - 1].append(n_targets + i)
                cols = [set(range(n_targets))]
                for i in range(min(nt, len(new))): cols.append({n_targets + i})
                for kd in sorted(set(kinds[i] for i in range(nt, len(new)))):
                    cols.append(set(n_targets + i for i in range(nt, len(new)) if kinds[i] == kd))
                g = pynauty.Graph(N, directed=False, vertex_coloring=cols); g.set_adjacency_dict(adj)
                ce = pynauty.certificate(g)
                if ce not in nxt: nxt[ce] = new
        states = nxt
    return list(states.values())
def bipartite_graphs(pdeg, gdeg):
    P = sorted(pdeg); G = sorted(gdeg); cap = dict(gdeg); out = []
    def rec(k, edges):
        if k == len(P):
            if all(v == 0 for v in cap.values()): out.append(tuple(edges))
            return
        j = P[k]; avail = [a for a in G if cap[a] > 0]
        for S in itertools.combinations(avail, pdeg[j]):
            for a in S: cap[a] -= 1
            rec(k + 1, edges + [(j, a) for a in S])
            for a in S: cap[a] += 1
    rec(0, [])
    return out
def job_cert(lay, pat, lat, edges):
    c, cp = lay.c, lay.cp
    vid = {}
    def v(x):
        if x not in vid: vid[x] = len(vid)
        return vid[x]
    for j in range(1, c + 1): v(('P', j))
    for a in range(1, cp + 1): v(('G', a))
    for k in range(len(pat)): v(('pair', k))
    for k in range(len(lat)): v(('pline', k))
    adj = {i: [] for i in range(len(vid))}
    def link(x, y): adj[v(x)].append(v(y)); adj[v(y)].append(v(x))
    for (j, a) in edges: link(('P', j), ('G', a))
    for k, pr in enumerate(pat):
        for a in pr: link(('pair', k), ('G', a))
    for k, pr in enumerate(lat):
        for j in pr: link(('pline', k), ('P', j))
    for k in range(lay.nt): link(('pair', k), ('pline', k))
    colours = [set(v(('P', j)) for j in range(1, c + 1)), set(v(('G', a)) for a in range(1, cp + 1))]
    for kind in sorted(set(('pair', k < lay.nt, lay.pair_axes[k]) for k in range(len(pat)))):
        colours.append(set(v(('pair', k)) for k in range(len(pat)) if (k < lay.nt, lay.pair_axes[k]) == kind[1:]))
    for kind in sorted(set(('pline', k < lay.nt, lay.line_centers[k]) for k in range(len(lat)))):
        colours.append(set(v(('pline', k)) for k in range(len(lat)) if (k < lay.nt, lay.line_centers[k]) == kind[1:]))
    g = pynauty.Graph(len(vid), directed=False, vertex_coloring=[cc for cc in colours if cc]); g.set_adjacency_dict(adj)
    return pynauty.certificate(g)
def jobs_for_case(lay):
    c, cp = lay.c, lay.cp
    pat_opts = attachments(cp, lay.pair_axes, [1] * lay.nt + [2] * lay.nb, lay.nt)
    lat_opts = attachments(c, lay.line_centers, [1] * lay.nt + [2] * lay.nbl, lay.nt)
    seen, jobs = set(), []
    for pat in pat_opts:
        gdeg = {a: 4 - sum(1 for pr in pat if a in pr) for a in range(1, cp + 1)}
        if any(d < 0 for d in gdeg.values()): continue
        for lat in lat_opts:
            pdeg = {j: 4 - sum(1 for pr in lat if j in pr) for j in range(1, c + 1)}
            if any(d < 0 for d in pdeg.values()) or sum(gdeg.values()) != sum(pdeg.values()): continue
            for edges in bipartite_graphs(pdeg, gdeg):
                ce = job_cert(lay, pat, lat, edges)
                if ce in seen: continue
                seen.add(ce); jobs.append((pat, lat, edges))
    return jobs
def enumerate_voltages(lay, pat, lat, edges):
    n, c, cp = lay.n, lay.c, lay.cp
    resP = {j: set(V) for j in range(1, c + 1)}; resG = {a: set(V) for a in range(1, cp + 1)}
    pat_cos, lat_cos = [], []
    for k, pr in enumerate(pat):
        ax = lay.pair_axes[k]; opts = []
        for y in pr:
            if any(al not in SUB[ax] for al in resG[y]): opts.append([E]); resG[y] = resG[y] & set(SUB[ax])
            else: opts.append([E, OTHER[ax]])
        pat_cos.append(opts)
    for k, pr in enumerate(lat):
        ce = lay.line_centers[k]; opts = []
        for y in pr:
            if any(al not in SUB[ce] for al in resP[y]): opts.append([E]); resP[y] = resP[y] & set(SUB[ce])
            else: opts.append([E, OTHER[ce]])
        lat_cos.append(opts)
    adj = {('P', j): [] for j in range(1, c + 1)}; adj.update({('G', a): [] for a in range(1, cp + 1)})
    for (j, a) in edges: adj[('P', j)].append(('G', a)); adj[('G', a)].append(('P', j))
    tree_child, visited, order = {}, set(), []
    for root in list(adj):
        if root in visited: continue
        visited.add(root); stack = [root]
        while stack:
            u = stack.pop()
            for w in adj[u]:
                if w not in visited:
                    visited.add(w); stack.append(w)
                    e = (u[1], w[1]) if u[0] == 'P' else (w[1], u[1]); tree_child[e] = w; order.append(e)
    edge_order = order + [e for e in edges if e not in tree_child]
    volt_opts = []
    for e in edge_order:
        if e in tree_child:
            w = tree_child[e]; H = resP[w[1]] if w[0] == 'P' else resG[w[1]]
            volt_opts.append(coset_reps(H))
            if w[0] == 'P': resP[w[1]] = {E}
            else: resG[w[1]] = {E}
        else: volt_opts.append(list(V))
    base = [0] * n
    def inc(p, l): base[p] |= 1 << l
    for i in lay.centers:
        for j in lay.centers:
            if i != j: inc(lay.O(i), lay.axis(j))
        for k, ce in enumerate(lay.tpairs):
            if ce == i:
                inc(lay.O(i), lay.pline(k, 0)); inc(lay.O(i), lay.pline(k, 1))
    for k, ax in enumerate(lay.tpairs):
        inc(lay.pair(k, 0), lay.axis(ax)); inc(lay.pair(k, 1), lay.axis(ax))
        inc(lay.pair(k, 0), lay.pline(k, 0)); inc(lay.pair(k, 1), lay.pline(k, 1))
    out = []
    for pc in itertools.product(*[itertools.product(*o) for o in pat_cos]):
        for lc in itertools.product(*[itertools.product(*o) for o in lat_cos]):
            pl = list(base)
            for k, pr in enumerate(pat):
                ax = lay.pair_axes[k]; sig = SIG[ax]; other = OTHER[ax]
                for (g, co) in zip(pr, pc[k]):
                    for al in (E, sig):
                        pl[lay.pair(k, 0)] |= 1 << lay.G(g, vmul(co, al))
                        pl[lay.pair(k, 1)] |= 1 << lay.G(g, vmul(vmul(co, al), other))
            for k, pr in enumerate(lat):
                ce = lay.line_centers[k]; sig = SIG[ce]; other = OTHER[ce]
                for (pj, co) in zip(pr, lc[k]):
                    for al in (E, sig):
                        pl[lay.P(pj, vmul(co, al))] |= 1 << lay.pline(k, 0)
                        pl[lay.P(pj, vmul(vmul(co, al), other))] |= 1 << lay.pline(k, 1)
            if not c4_ok_all(pl, n): continue
            def rec(k, cur):
                if k == len(edge_order):
                    if all(bin(m).count('1') == 4 for m in cur):
                        lp = [0] * n
                        for p in range(n):
                            m = cur[p]
                            while m:
                                l = (m & -m).bit_length() - 1; lp[l] |= 1 << p; m &= m - 1
                        if all(bin(m).count('1') == 4 for m in lp): out.append(tuple(cur))
                    return
                (j, a) = edge_order[k]
                for v in volt_opts[k]:
                    nxt = list(cur); ok = True; touched = []
                    for be in V:
                        p = lay.P(j, be); nxt[p] |= 1 << lay.G(a, vmul(be, v)); touched.append(p)
                    for p in touched:
                        mp = nxt[p]
                        for q in range(n):
                            if q != p and bin(mp & nxt[q]).count('1') > 1: ok = False; break
                        if not ok: break
                    if ok: rec(k + 1, nxt)
            rec(0, pl)
    return out
def c4_ok_all(pl, n):
    for p in range(n):
        for q in range(p + 1, n):
            if bin(pl[p] & pl[q]).count('1') > 1: return False
    return True
def certificate(lay, pl):
    n = lay.n; N0 = 2 * n; adjd = {i: [] for i in range(N0)}
    for p in range(n):
        m = pl[p]
        while m:
            l = (m & -m).bit_length() - 1; adjd[p].append(n + l); adjd[n + l].append(p); m &= m - 1
    aux = []
    for si, s in ((1, S1), (2, S2), (3, S3)):
        for p in range(n):
            q = lay.act_pt(s, p)
            if q > p: aux.append((si, p, q))
        for l in range(n):
            k = lay.act_ln(s, l)
            if k > l: aux.append((si, n + l, n + k))
    for idx, (si, u, w) in enumerate(aux): adjd[N0 + idx] = [u, w]
    colours = [set(range(n)), set(range(n, N0))] + [set(N0 + i for i, (si, u, w) in enumerate(aux) if si == s) for s in (1, 2, 3)]
    g = pynauty.Graph(N0 + len(aux), directed=False, vertex_coloring=[cc for cc in colours if cc]); g.set_adjacency_dict(adjd)
    return pynauty.certificate(g)
def pt_on_axis(i, var, sgn):
    return {1: ((0, None), (1, None), (sgn, var)), 2: ((1, None), (0, None), (sgn, var)), 3: ((1, None), (sgn, var), (0, None))}[i]
def line_through_centre(i, var, sgn):
    return {1: ((0, None), (1, None), (sgn, var)), 2: ((1, None), (0, None), (sgn, var)), 3: ((1, None), (sgn, var), (0, None))}[i]
def tline(i, var):
    return {1: (((0, None), (1, var), (-1, None)), ((0, None), (1, var), (1, None))),
            2: (((1, var), (0, None), (-1, None)), ((1, var), (0, None), (1, None))),
            3: (((1, var), (-1, None), (0, None)), ((1, var), (1, None), (0, None)))}[i]
def sym_coords(lay):
    pts, lns = {}, {}
    unit = lambda i: tuple((1, None) if k == i - 1 else (0, None) for k in range(3))
    for i in lay.centers: pts[lay.O(i)] = unit(i); lns[lay.axis(i)] = unit(i)
    for k, ax in enumerate(lay.tpairs):
        pts[lay.pair(k, 0)] = pt_on_axis(ax, 't%d' % k, 1); pts[lay.pair(k, 1)] = pt_on_axis(ax, 't%d' % k, -1)
        lns[lay.pline(k, 0)], lns[lay.pline(k, 1)] = tline(ax, 't%d' % k)
    for k, ax in enumerate(lay.bpairs):
        kk = lay.nt + k; pts[lay.pair(kk, 0)] = pt_on_axis(ax, 'S%d' % k, 1); pts[lay.pair(kk, 1)] = pt_on_axis(ax, 'S%d' % k, -1)
    for k, ce in enumerate(lay.bplines):
        kk = lay.nt + k; lns[lay.pline(kk, 0)] = line_through_centre(ce, 'T%d' % k, 1); lns[lay.pline(kk, 1)] = line_through_centre(ce, 'T%d' % k, -1)
    for j in range(1, lay.c + 1):
        for al in V:
            u = (al[0] * al[1], 'u%d' % j if j > 1 else None); v = (al[0] * al[2], 'v%d' % j if j > 1 else None)
            pts[lay.P(j, al)] = ((1, None), u, v)
    for a in range(1, lay.cp + 1):
        for al in V: lns[lay.G(a, al)] = ((1, None), (al[0] * al[1], 'm%d' % a), (al[0] * al[2], 'n%d' % a))
    return pts, lns
def variables(lay):
    return (['t%d' % k for k in range(lay.nt)] + ['S%d' % k for k in range(lay.nb)] + ['u%d' % j for j in range(2, lay.c + 1)]
            + ['v%d' % j for j in range(2, lay.c + 1)] + ['T%d' % k for k in range(lay.nbl)]
            + ['m%d' % a for a in range(1, lay.cp + 1)] + ['n%d' % a for a in range(1, lay.cp + 1)])
def dotpoly(P, L):
    terms = {}
    for (s1, x1), (s2, x2) in zip(P, L):
        if s1 == 0 or s2 == 0: continue
        mono = tuple(sorted(v for v in (x1, x2) if v is not None)); terms[mono] = terms.get(mono, 0) + s1 * s2
    terms = {m: c for m, c in terms.items() if c != 0}
    if not terms: return None
    items = sorted(terms.items())
    if items[0][1] < 0: items = [(m, -c) for m, c in items]
    return tuple(items)
def poly_str(p):
    return '+'.join(('%d' % c) if not mono else (('' if c == 1 else '-' if c == -1 else '%d*' % c) + '*'.join(mono)) for mono, c in p).replace('+-', '-')
def build(lay, pl):
    pts, lns = sym_coords(lay); eqs, ineqs = set(), set(); n = lay.n
    for p in range(n):
        for l in range(n):
            poly = dotpoly(pts[p], lns[l])
            if pl[p] >> l & 1:
                if poly is None: continue
                if len(poly) == 1 and poly[0][0] == (): return None, None
                eqs.add(poly)
            else:
                if poly is None: return None, None
                if len(poly) == 1 and poly[0][0] == (): continue
                ineqs.add(poly)
    return sorted(eqs), sorted(ineqs)
def singular_script(lay, eqs, tag):
    vs = variables(lay)
    s = 'ring r = 0,(%s),dp;\nideal I = %s;\nideal J = std(I);\nint d0 = dim(J);\nif (d0 >= 0) {\n' % (','.join(vs), ','.join(poly_str(e) for e in eqs))
    for v in vs: s += '  J = sat(J, ideal(%s));\n' % v
    pvars = [('t%d' % k, ax) for k, ax in enumerate(lay.tpairs)] + [('S%d' % k, ax) for k, ax in enumerate(lay.bpairs)]
    for (x, ax), (y, ay) in itertools.combinations(pvars, 2):
        if ax == ay: s += '  J = sat(J, ideal(%s^2-%s^2));\n' % (x, y)
    lvars = [('t%d' % k, ce) for k, ce in enumerate(lay.tpairs)] + [('T%d' % k, ce) for k, ce in enumerate(lay.bplines)]
    for (x, cx), (y, cy) in itertools.combinations(lvars, 2):
        if cx == cy and not (x.startswith('t') and y.startswith('t')): s += '  J = sat(J, ideal(%s^2-%s^2));\n' % (x, y)
    for a, b in itertools.combinations(range(1, lay.c + 1), 2):
        ua = 'u%d' % a if a > 1 else '1'; va = 'v%d' % a if a > 1 else '1'
        s += '  J = sat(J, ideal(%s^2-u%d^2, %s^2-v%d^2));\n' % (ua, b, va, b)
    for a, b in itertools.combinations(range(1, lay.cp + 1), 2): s += '  J = sat(J, ideal(m%d^2-m%d^2, n%d^2-n%d^2));\n' % (a, b, a, b)
    s += '  J = std(J);\n}\nint d = dim(J);\n'
    s += 'print("RES %s dim0=" + string(d0) + " dim=" + string(d) + " vdim=" + string(vdim(J)));\n' % tag
    s += 'if (d == 0 and vdim(J) > 0) { def R = solve(J, 40, 0, "nodisplay"); setring R; print("SOLS %s"); print(SOL); print("ENDSOLS"); kill R; }\nkill r;\n' % tag
    return s
def parse_complex(s):
    s = s.strip().replace(' ', ''); s2 = re.sub(r'(?<![0-9a-zA-Z_.])i(?![0-9a-zA-Z_])', '1j', s)
    return complex(eval(s2, {'__builtins__': {}}, {}))
def parse_output(out, vs):
    res = {}
    for m in re.finditer(r'RES (\S+) dim0=(-?\d+) dim=(-?\d+) vdim=(-?\d+)', out):
        res[m.group(1)] = {'dim0': int(m.group(2)), 'dim': int(m.group(3)), 'vdim': int(m.group(4)), 'sols': []}
    for m in re.finditer(r'SOLS (\S+)\n(.*?)ENDSOLS', out, flags=re.S):
        vals = [l.strip() for l in m.group(2).splitlines() if l.strip() and not re.match(r'^\[\d+\]:$', l.strip())]
        res[m.group(1)]['sols'] = [{vs[j]: parse_complex(vals[k + j]) for j in range(len(vs))} for k in range(0, len(vals) - len(vs) + 1, len(vs))]
    return res
def is_real(sol): return all(abs(v.imag) < 1e-12 * max(1.0, abs(v)) for v in sol.values())
def cases(n):
    out = []
    def vecs(total, cap):
        return [b for b in itertools.product(range(cap + 1), repeat=3) if sum(b) == total]
    if n % 2 == 0:
        for c in range(1, n // 4 + 1):
            for cp in range(1, n // 4 + 1):
                sb, sbp = n // 2 - 2 * c, n // 2 - 2 * cp
                if sb < 0 or sbp < 0: continue
                seen = set()
                for b in vecs(sb, cp // 2):
                    for bp in vecs(sbp, c // 2):
                        key = min((tuple(b[i] for i in perm), tuple(bp[i] for i in perm)) for perm in itertools.permutations(range(3)))
                        if key in seen: continue
                        seen.add(key)
                        bl = [i + 1 for i in range(3) for _ in range(key[0][i])]; bpl = [i + 1 for i in range(3) for _ in range(key[1][i])]
                        out.append(('A', c, cp, bl, bpl))
    else:
        for c in range(2, (n - 5) // 4 + 1):
            for cp in range(2, (n - 5) // 4 + 1):
                sb, sbp = (n - 5) // 2 - 2 * c, (n - 5) // 2 - 2 * cp
                if sb < 0 or sbp < 0: continue
                seen = set()
                for b2, b3 in [(x, sb - x) for x in range(sb + 1)]:
                    if max(b2, b3) > cp // 2: continue
                    for bp2, bp3 in [(x, sbp - x) for x in range(sbp + 1)]:
                        if max(bp2, bp3) > c // 2: continue
                        key = min(((b2, b3), (bp2, bp3)), ((b3, b2), (bp3, bp2)))
                        if key in seen: continue
                        seen.add(key)
                        bl = [2] * key[0][0] + [3] * key[0][1]; bpl = [2] * key[1][0] + [3] * key[1][1]
                        out.append(('B', c, cp, bl, bpl))
        if n % 4 == 1 and n >= 13:
            out.append(('C', (n - 9) // 4, (n - 9) // 4, [], []))
    return out
def tag_of(n, cs): return '%d_%s_%d_%d_%s_%s' % (n, cs[0], cs[1], cs[2], ''.join(map(str, cs[3])) or '-', ''.join(map(str, cs[4])) or '-')
if __name__ == '__main__':
    mode = sys.argv[1]; n = int(sys.argv[2]) if mode != 'solve' else None
    if mode == 'cases':
        for cs in cases(n): print(tag_of(n, cs))
    elif mode == 'enum':
        idx = int(sys.argv[3]); cs = cases(n)[idx]; lay = Layout(n, cs[0], cs[1], cs[2], cs[3], cs[4]); t0 = time.time()
        fn = '/home/claude/kn/%s.json' % tag_of(n, cs)
        jobs = jobs_for_case(lay); tj = time.time() - t0
        try: data = json.load(open(fn)); seen = {bytes.fromhex(k): v for k, v in data['structures'].items()}; done = data.get('jobs_done', 0)
        except Exception: seen, done = {}, 0
        start = done; last = start
        for jidx in range(start, len(jobs)):
            pat, lat, edges = jobs[jidx]
            for pl in enumerate_voltages(lay, pat, lat, edges):
                ce = certificate(lay, pl)
                if ce not in seen: seen[ce] = pl
            last = jidx + 1
            if time.time() - t0 > 80: break
        json.dump({'n': n, 'case': cs, 'structures': {k.hex(): v for k, v in seen.items()}, 'complete': last >= len(jobs), 'njobs': len(jobs), 'jobs_done': last}, open(fn, 'w'))
        print("%s: jobs %d (%.1fs to list), done %d, structures %d (%.1fs)%s" % (tag_of(n, cs), len(jobs), tj, last, len(seen), time.time() - t0, '' if last >= len(jobs) else '  PARTIAL'), flush=True)
    elif mode == 'solve':
        fn = sys.argv[2]; data = json.load(open(fn)); cs = data['case']; n = data['n']; lay = Layout(n, cs[0], cs[1], cs[2], cs[3], cs[4]); vs = variables(lay)
        keys = sorted(data['structures']); script = 'LIB "solve.lib";\n'; tags = []
        rf = fn.replace('.json', '_results.json')
        try: results = json.load(open(rf))
        except Exception: results = {}
        CH = int(sys.argv[3]) if len(sys.argv) > 3 else 4000; todo = [k for k in keys if k not in results][:CH]
        for idx, key in enumerate(keys):
            if key not in todo: continue
            eqs, ineqs = build(lay, data['structures'][key])
            if eqs is None: results[key] = {'status': 'contradiction'}; continue
            script += singular_script(lay, eqs, 'S%d' % idx); tags.append((idx, key, ineqs))
        script += 'quit;\n'; open('/tmp/kn.sing', 'w').write(script); t0 = time.time()
        try: out = subprocess.run(['Singular', '-q', '/tmp/kn.sing'], capture_output=True, text=True, timeout=90).stdout
        except subprocess.TimeoutExpired as ex: out = (ex.stdout or b'').decode() if isinstance(ex.stdout, bytes) else (ex.stdout or ''); print("SINGULAR TIMEOUT", fn)
        parsed = parse_output(out, vs)
        for idx, key, ineqs in tags:
            r = parsed.get('S%d' % idx)
            if r is None: results[key] = {'status': 'missing'}; continue
            good = []
            for sol in r['sols']:
                if any(abs(eval(poly_str(q).replace('^', '**'), {}, dict(sol))) < 1e-12 for q in ineqs): continue
                good.append({'real': is_real(sol), 'sol': {k: [v.real, v.imag] for k, v in sol.items()}})
            results[key] = {'status': 'ok', 'dim': r['dim'], 'vdim': r['vdim'], 'nsols': len(r['sols']), 'nondeg': len(good), 'nreal': sum(1 for g in good if g['real']), 'sols': good}
        json.dump(results, open(rf, 'w'))
        if len([k for k in keys if k not in results]) > 0: print('  (%d structures still unsolved)' % len([k for k in keys if k not in results]))
        ne = [(i, results[k]['dim'], results[k]['vdim'], results[k]['nondeg'], results[k]['nreal']) for i, k in enumerate(keys) if results[k].get('status') == 'ok' and results[k]['dim'] >= 0]
        print("%s: %d structures, %.1fs; nonempty: %s" % (fn.split('/')[-1], len(keys), time.time() - t0, ne if len(ne) < 12 else (len(ne), 'structures;', sum(x[3] for x in ne), 'nondeg,', sum(x[4] for x in ne), 'real')), flush=True)
