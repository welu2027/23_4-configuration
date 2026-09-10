import itertools, json, sys, time, re, subprocess
import pynauty
V = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
E, S1, S2, S3 = V
VIDX = {v: i for i, v in enumerate(V)}
def vmul(a, b): return tuple(x * y for x, y in zip(a, b))
SUB = {1: (E, S1), 2: (E, S2), 3: (E, S3)}
SIG = {1: S1, 2: S2, 3: S3}
class Layout:
    def __init__(self, c, cprime, b_axis, bp_centers):
        self.c, self.cp, self.b_axis, self.bp = c, cprime, b_axis, bp_centers
        self.np = 7 + 4 * c; assert self.np == 23
        self.gstart = 5 + 2 * len(bp_centers)
        self.nl = self.gstart + 4 * cprime; assert self.nl == 23
    def T(self, i, pr): return 1 + 2 * (i - 1) + pr
    def B(self, pr): return 5 + pr
    def P(self, j, a): return 7 + 4 * (j - 1) + VIDX[a]
    def Tp(self, i, pr): return 1 + 2 * (i - 1) + pr
    def Bp(self, k, pr): return 5 + 2 * k + pr
    def G(self, a, al): return self.gstart + 4 * (a - 1) + VIDX[al]
    def act_pt(self, s, v):
        if v == 0: return 0
        if 1 <= v <= 4: return v if s in (E, S1) else (v + 1 if v % 2 == 1 else v - 1)
        if v in (5, 6): return v if s in (E, SIG[self.b_axis]) else (11 - v)
        j = (v - 7) // 4 + 1; al = V[(v - 7) % 4]; return self.P(j, vmul(al, s))
    def act_ln(self, s, v):
        if v == 0: return 0
        if 1 <= v <= 4: return v if s in (E, S1) else (v + 1 if v % 2 == 1 else v - 1)
        if v < self.gstart:
            k = (v - 5) // 2; pr = (v - 5) % 2; cen = self.bp[k]
            return v if s in (E, SIG[cen]) else self.Bp(k, 1 - pr)
        a = (v - self.gstart) // 4 + 1; al = V[(v - self.gstart) % 4]; return self.G(a, vmul(al, s))
def coset_reps(H):
    reps, seen = [], set()
    for v in V:
        if v in seen: continue
        reps.append(v); seen |= {vmul(v, h) for h in H}
    return reps
def canon_key(lay, bxy, bpch):
    c, cp = lay.c, lay.cp
    permsP = [{1:1,2:2,3:3,4:4}, {1:2,2:1,3:3,4:4}, {1:1,2:2,3:4,4:3}, {1:2,2:1,3:4,4:3}]
    keys = []
    for pP in permsP:
        swap12 = (pP[1] == 2)
        gperms = [{a: a for a in range(1, cp + 1)}]
        if cp == 4: gperms.append({1: 1, 2: 2, 3: 4, 4: 3})
        for pG in gperms:
            pG2 = dict(pG)
            if swap12: pG2[1], pG2[2] = pG[2], pG[1]
            bx = tuple(sorted(pG2[a] for a in bxy))
            imgs = [tuple(sorted(pP[j] for j in pr)) for pr in bpch]
            grouped = []
            for cen in sorted(set(lay.bp)):
                grouped.append(tuple(sorted(imgs[k] for k in range(len(lay.bp)) if lay.bp[k] == cen)))
            keys.append((bx, tuple(grouped)))
    return min(keys)
def jobs_for_case(lay):
    c, cp, nB = lay.c, lay.cp, len(lay.bp)
    P_orbs = list(range(1, c + 1)); G_orbs = list(range(1, cp + 1))
    pair_lists = [list(itertools.combinations(P_orbs, 2)) for _ in range(nB)]
    seen, jobs = set(), []
    for bxy in itertools.combinations(G_orbs, 2):
        for bpch in itertools.product(*pair_lists):
            ok = True
            for j in P_orbs:
                for cen in set(lay.bp):
                    if sum(1 for k in range(nB) if lay.bp[k] == cen and j in bpch[k]) > 1: ok = False
            if not ok: continue
            key = canon_key(lay, bxy, bpch)
            if key in seen: continue
            seen.add(key)
            gdeg = {a: 4 - (1 if a <= 2 else 0) - (1 if a in bxy else 0) for a in G_orbs}
            pdeg = {j: 4 - (1 if j <= 2 else 0) - sum(1 for k in range(nB) if j in bpch[k]) for j in P_orbs}
            if any(d < 0 or d > cp for d in pdeg.values()) or any(d < 0 or d > c for d in gdeg.values()): continue
            if sum(pdeg.values()) != sum(gdeg.values()): continue
            Eg = sum(pdeg.values())
            allpairs = [(j, a) for j in P_orbs for a in G_orbs]
            for edges in itertools.combinations(allpairs, Eg):
                if any(sum(1 for (j, a) in edges if j == jj) != pdeg[jj] for jj in P_orbs): continue
                if any(sum(1 for (j, a) in edges if a == aa) != gdeg[aa] for aa in G_orbs): continue
                jobs.append((bxy, bpch, edges))
    return jobs
def enumerate_voltages(lay, bxy, bpch, edges):
    c, cp, nB = lay.c, lay.cp, len(lay.bp)
    resP = {j: set(V) for j in range(1, c + 1)}; resG = {a: set(V) for a in range(1, cp + 1)}
    for j in (1, 2): resP[j] = set(SUB[1])
    for a in (1, 2): resG[a] = set(SUB[1])
    pl = [0] * 23
    def inc(p, l):
        pl[p] |= 1 << l
    for pr in (0, 1):
        for i in (1, 2): inc(0, lay.Tp(i, pr)); inc(lay.T(i, pr), 0); inc(lay.T(i, pr), lay.Tp(i, pr))
    for i in (1, 2):
        for al in (E, S1): inc(lay.T(i, 0), lay.G(i, al))
        for al in (S2, S3): inc(lay.T(i, 1), lay.G(i, al))
    for i in (1, 2):
        for al in (E, S1): inc(lay.P(i, al), lay.Tp(i, 0))
        for al in (S2, S3): inc(lay.P(i, al), lay.Tp(i, 1))
    sB = SIG[lay.b_axis]
    def use_freedom(res, orb, H):
        if any(al not in H for al in res[orb]):
            res[orb] = res[orb] & set(H); return [E]
        return [E, S1]
    bcos = []
    x, y = bxy
    bcos.append(use_freedom(resG, x, SUB[lay.b_axis]))
    bcos.append(use_freedom(resG, y, SUB[lay.b_axis]))
    bpcos = []
    for k in range(nB):
        cen = lay.bp[k]; xp, yp = bpch[k]
        opts_k = [use_freedom(resP, xp, SUB[cen]), use_freedom(resP, yp, SUB[cen])]
        bpcos.append(opts_k)
    adj = {('P', j): [] for j in range(1, c + 1)}; adj.update({('G', a): [] for a in range(1, cp + 1)})
    for (j, a) in edges: adj[('P', j)].append(('G', a)); adj[('G', a)].append(('P', j))
    tree_child = {}
    visited = set(); order = []
    for root in list(adj):
        if root in visited: continue
        visited.add(root); stack = [root]
        while stack:
            u = stack.pop()
            for w in adj[u]:
                if w not in visited:
                    visited.add(w); stack.append(w)
                    e = (u[1], w[1]) if u[0] == 'P' else (w[1], u[1])
                    tree_child[e] = w; order.append(e)
    nontree = [e for e in edges if e not in tree_child]
    edge_order = order + nontree
    volt_opts = []
    for e in edge_order:
        if e in tree_child:
            w = tree_child[e]
            H = resP[w[1]] if w[0] == 'P' else resG[w[1]]
            volt_opts.append(coset_reps(H))
            if w[0] == 'P': resP[w[1]] = {E}
            else: resG[w[1]] = {E}
        else:
            volt_opts.append(list(V))
    out = []
    lp_cache = None
    for bc in itertools.product(*bcos):
        for bpc in itertools.product(*[itertools.product(*o) for o in bpcos]):
            pl0 = list(pl)
            for (g, co) in zip(bxy, bc):
                for al in (E, sB):
                    pl0[lay.B(0)] |= 1 << lay.G(g, vmul(co, al))
                    pl0[lay.B(1)] |= 1 << lay.G(g, vmul(vmul(co, al), S1))
            for k in range(nB):
                cen = lay.bp[k]
                for (pj, co) in zip(bpch[k], bpc[k]):
                    for al in (E, SIG[cen]):
                        pl0[lay.P(pj, vmul(co, al))] |= 1 << lay.Bp(k, 0)
                        pl0[lay.P(pj, vmul(vmul(co, al), S1))] |= 1 << lay.Bp(k, 1)
            if not c4_ok_all(pl0): continue
            def rec(k, cur):
                if k == len(edge_order):
                    if all(bin(m).count('1') == 4 for m in cur):
                        lp = [0] * 23
                        for p in range(23):
                            m = cur[p]
                            while m:
                                l = (m & -m).bit_length() - 1; lp[l] |= 1 << p; m &= m - 1
                        if all(bin(m).count('1') == 4 for m in lp):
                            out.append(tuple(cur))
                    return
                (j, a) = edge_order[k]
                for v in volt_opts[k]:
                    nxt = list(cur); ok = True; touched = []
                    for be in V:
                        p = lay.P(j, be); l = lay.G(a, vmul(be, v))
                        nxt[p] |= 1 << l; touched.append(p)
                    for p in touched:
                        mp = nxt[p]
                        for q in range(23):
                            if q != p and bin(mp & nxt[q]).count('1') > 1: ok = False; break
                        if not ok: break
                    if ok: rec(k + 1, nxt)
            rec(0, pl0)
    return out
def c4_ok_all(pl):
    for p in range(23):
        for q in range(p + 1, 23):
            if bin(pl[p] & pl[q]).count('1') > 1: return False
    return True
def certificate(lay, pl):
    n = 46
    adjd = {i: [] for i in range(n)}
    for p in range(23):
        m = pl[p]
        while m:
            l = (m & -m).bit_length() - 1; adjd[p].append(23 + l); adjd[23 + l].append(p); m &= m - 1
    aux = []
    for si, s in ((1, S1), (2, S2), (3, S3)):
        for p in range(23):
            q = lay.act_pt(s, p)
            if q > p: aux.append((si, p, q))
        for l in range(23):
            k = lay.act_ln(s, l)
            if k > l: aux.append((si, 23 + l, 23 + k))
    N = n + len(aux)
    for idx, (si, u, w) in enumerate(aux): adjd[n + idx] = [u, w]
    colours = [set(range(23)), set(range(23, 46))] + [set(n + i for i, (si, u, w) in enumerate(aux) if si == s) for s in (1, 2, 3)]
    g = pynauty.Graph(N, directed=False, vertex_coloring=[cc for cc in colours if cc])
    g.set_adjacency_dict(adjd)
    return pynauty.certificate(g)
if __name__ == '__main__':
    cp = int(sys.argv[1]); b_axis = int(sys.argv[2]); bp = [int(x) for x in sys.argv[3].split(',')]
    lay = Layout(4, cp, b_axis, bp)
    jobs = jobs_for_case(lay)
    tag = '%d_%d_%s' % (cp, b_axis, ''.join(map(str, bp)))
    if len(sys.argv) == 4:
        print("case", tag, "jobs:", len(jobs)); sys.exit(0)
    start, end = int(sys.argv[4]), int(sys.argv[5])
    t0 = time.time(); seen = {}
    fn = '/home/claude/genVfix_%s_%d_%d.json' % (tag, start, end)
    for idx in range(start, min(end, len(jobs))):
        bxy, bpch, edges = jobs[idx]
        res = enumerate_voltages(lay, bxy, bpch, edges)
        for pl in res:
            ce = certificate(lay, pl)
            if ce not in seen: seen[ce] = pl
        if time.time() - t0 > 80:
            print("TIME LIMIT at job", idx); end = idx + 1; break
    json.dump({'cp': cp, 'b_axis': b_axis, 'bp': bp, 'jobs_done': [start, end], 'njobs': len(jobs),
               'structures': {ce.hex(): pl for ce, pl in seen.items()}}, open(fn, 'w'))
    print("case %s jobs %d-%d of %d: distinct structures %d (%.1fs)" % (tag, start, end, len(jobs), len(seen), time.time() - t0))
