import json, sys, time, pynauty
from classify_Vn import *
def coloured_graph(lay, pl):
    n = lay.n; N0 = 2 * n; adj = {i: set() for i in range(N0)}
    for p in range(n):
        m = pl[p]
        while m:
            l = (m & -m).bit_length() - 1; adj[p].add(n + l); adj[n + l].add(p); m &= m - 1
    aux = []
    for si, s in ((1, S1), (2, S2), (3, S3)):
        for p in range(n):
            q = lay.act_pt(s, p)
            if q > p: aux.append((si, p, q))
        for l in range(n):
            k = lay.act_ln(s, l)
            if k > l: aux.append((si, n + l, n + k))
    for idx, (si, u, w) in enumerate(aux): adj[N0 + idx] = {u, w}; adj[u].add(N0 + idx); adj[w].add(N0 + idx)
    colour = [0] * n + [1] * n + [1 + si for (si, u, w) in aux]
    return adj, colour
def wl_hash(adj, colour, rounds=6):
    col = list(colour)
    for _ in range(rounds):
        sig = [(col[v], tuple(sorted(col[w] for w in adj[v]))) for v in range(len(col))]
        table = {s: i for i, s in enumerate(sorted(set(sig)))}; col = [table[s] for s in sig]
    return hash(tuple(sorted(col)))
def relabelled(adj, colour, lab):
    edges = frozenset(frozenset((lab[u], lab[w])) for u in adj for w in adj[u] if u < w)
    cols = tuple(sorted((lab[v], colour[v]) for v in range(len(colour))))
    return edges, cols
case = sys.argv[1]; data = json.load(open('kn/%s.json' % case)); cs = data['case']; lay = Layout(23, cs[0], cs[1], cs[2], cs[3], cs[4])
t0 = time.time(); jobs = jobs_for_case(lay); raw = []
for (pat, lat, edges) in jobs:
    raw.extend(enumerate_voltages(lay, pat, lat, edges))
print("%s: raw structures %d from %d jobs (%.1fs)" % (case, len(raw), len(jobs), time.time() - t0), flush=True)
classes = {}; graphs = []
for i, pl in enumerate(raw):
    adj, colour = coloured_graph(lay, pl); graphs.append((adj, colour))
    classes.setdefault(certificate(lay, pl), []).append(i)
print("nauty classes:", len(classes), " (paper's count for this case: %d)" % len(data['structures']))
wl = [wl_hash(*g) for g in graphs]
wl_classes = {}
for i, h in enumerate(wl): wl_classes.setdefault(h, []).append(i)
split = sum(1 for members in classes.values() if len({wl[i] for i in members}) > 1)
print("distinct WL hashes: %d ; nauty classes whose members have different WL hashes (would prove an over-merge): %d" % (len(wl_classes), split))
bad = 0; checked = 0
for cert, members in classes.items():
    if len(members) == 1: continue
    ref = members[0]; adj0, col0 = graphs[ref]
    g0 = pynauty.Graph(len(col0), directed=False, vertex_coloring=[set(v for v in range(len(col0)) if col0[v] == c) for c in sorted(set(col0))]); g0.set_adjacency_dict({v: sorted(adj0[v]) for v in adj0})
    can0 = relabelled(adj0, col0, {v: i for i, v in enumerate(pynauty.canon_label(g0))})
    for other in members[1:]:
        adj1, col1 = graphs[other]
        g1 = pynauty.Graph(len(col1), directed=False, vertex_coloring=[set(v for v in range(len(col1)) if col1[v] == c) for c in sorted(set(col1))]); g1.set_adjacency_dict({v: sorted(adj1[v]) for v in adj1})
        can1 = relabelled(adj1, col1, {v: i for i, v in enumerate(pynauty.canon_label(g1))})
        checked += 1
        if can0 != can1: bad += 1
print("merges checked with an explicit witness isomorphism: %d ; failures: %d  (%.1fs)" % (checked, bad, time.time() - t0))
open('/home/claude/dedup_audit_log.txt', 'a').write("%s raw=%d nauty=%d paper=%d wl=%d split=%d merges=%d bad=%d\n" % (case, len(raw), len(classes), len(data['structures']), len(wl_classes), split, checked, bad))
