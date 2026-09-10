from fractions import Fraction as F
from itertools import combinations, product
import json, pynauty

D=3
def mk(a,b=0): return (F(a),F(b))
def add(x,y): return (x[0]+y[0],x[1]+y[1])
def sub(x,y): return (x[0]-y[0],x[1]-y[1])
def mul(x,y): return (x[0]*y[0]+D*x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def iszero(x): return x[0]==0 and x[1]==0
r3=mk(0,1)
def q(a,b): return (F(a),F(b))

P0=[(q(0,0),q(1,0),q(2,1)), (q(1,0),q(0,0),q(1,F(-1,3))), (q(1,0),q(-1,F(1,3)),q(0,0)),
    (q(1,0),q(1,0),q(1,0)), (q(1,0),q(2,-1),q(2,-1)), (q(1,0),q(-1,F(2,3)),q(1,0)), (q(1,0),q(2,-1),q(0,F(-1,3)))]
L0=[(q(0,0),q(1,0),q(-1,0)), (q(1,0),q(0,0),q(-1,0)), (q(1,0),q(-2,-1),q(0,0)),
    (q(1,0),q(F(9,2),F(5,2)),q(F(-3,2),F(-1,2))), (q(1,0),q(F(3,2),F(1,2)),q(F(-3,2),F(1,2))),
    (q(1,0),q(F(1,2),F(1,2)),q(F(-3,2),F(-1,2))), (q(1,0),q(F(3,2),F(1,2)),q(F(-1,2),F(-1,2)))]
V=[(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]
def act(v,p): return tuple(mul(mk(s),c) for s,c in zip(v,p))
def cross(u,v):
    return (sub(mul(u[1],v[2]),mul(u[2],v[1])),sub(mul(u[2],v[0]),mul(u[0],v[2])),sub(mul(u[0],v[1]),mul(u[1],v[0])))
def same(u,v): return all(iszero(t) for t in cross(u,v))
def orbit(p):
    out=[]
    for v in V:
        r=act(v,p)
        if not any(same(r,s) for s in out): out.append(r)
    return out
P=[r for p in P0 for r in orbit(p)]; L=[r for l in L0 for r in orbit(l)]
print("points:",len(P)," lines:",len(L))
dot=lambda p,l: add(add(mul(p[0],l[0]),mul(p[1],l[1])),mul(p[2],l[2]))
inc=[[iszero(dot(p,l)) for l in L] for p in P]
ok=(all(sum(r)==4 for r in inc) and all(sum(inc[i][j] for i in range(22))==4 for j in range(22))
    and all(sum(a and b for a,b in zip(inc[i],inc[j]))<=1 for i,j in combinations(range(22),2))
    and all(not same(P[i],P[j]) for i,j in combinations(range(22),2))
    and all(not same(L[i],L[j]) for i,j in combinations(range(22),2)))
print("T3 is a geometric (22_4):",ok)
n=22
adj={i:[n+j for j in range(n) if inc[i][j]] for i in range(n)}
for j in range(n): adj[n+j]=[i for i in range(n) if inc[i][j]]
g=pynauty.Graph(2*n,adjacency_dict=adj,vertex_coloring=[set(range(n)),set(range(n,2*n))])
cert=pynauty.certificate(g).hex(); aut=pynauty.autgrp(g)[1]
g2=pynauty.Graph(2*n,adjacency_dict=adj); aut2=pynauty.autgrp(g2)[1]
print("T3 Levi aut kept apart:",aut," with dualities:",aut2," self-dual:",aut2>aut)
found={}
for fn in ['types_22_11_13_17.json','types_22_19.json']:
    for k,v in json.load(open(fn)).items():
        found.setdefault(k,[]).extend(v)
for k,occ in found.items():
    ps=sorted({o['p'] for o in occ})
    if k==cert: print("T3 MATCHES finite-field type with primes",ps,"aut",occ[0]['aut'],"/",occ[0]['aut_dual'])
multi=[(k,sorted({o['p'] for o in occ}),occ[0]['aut'],occ[0]['aut_dual']) for k,occ in found.items() if len({o['p'] for o in occ})>1]
print("multi-prime finite-field types:",[(k[:10],ps,a,b) for k,ps,a,b in multi])
print("T3 cert prefix:",cert[:10])
cinc=json.load(open('cuntz22_inc.json'))
adjc={i:[n+j for j in range(n) if cinc[i][j]] for i in range(n)}
for j in range(n): adjc[n+j]=[i for i in range(n) if cinc[i][j]]
gc=pynauty.Graph(2*n,adjacency_dict=adjc,vertex_coloring=[set(range(n)),set(range(n,2*n))])
cc=pynauty.certificate(gc).hex()
print("Cuntz cert prefix:",cc[:10]," Cuntz type among finite-field types:",cc in found, sorted({o['p'] for o in found.get(cc,[])}))
print("T3 == Cuntz type:",cert==cc)

import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
import sympy as sp
G=nx.Graph()
for i in range(22): G.add_node(('p',i),s=0); G.add_node(('l',i),s=1)
for i in range(22):
    for j in range(22):
        if inc[i][j]: G.add_edge(('p',i),('l',j))
swap=lambda a,b: a['s']!=b['s']
R3=sp.sqrt(3)
def tosp(x): return sp.Rational(x[0])+sp.Rational(x[1])*R3
Psp=[sp.Matrix([tosp(c) for c in p]) for p in P]; Lsp=[sp.Matrix([tosp(c) for c in l]) for l in L]
b=sp.symbols('b0:6'); Bm=sp.Matrix([[b[0],b[1],b[2]],[b[1],b[3],b[4]],[b[2],b[4],b[5]]])
found=0; tried=0
for m in GraphMatcher(G,G,node_match=swap).isomorphisms_iter():
    tried+=1
    eqs=[]
    for i in range(22):
        j=m[('p',i)][1]
        v=Bm*Psp[i]; w=Lsp[j]
        eqs+=[sp.expand(v[0]*w[1]-v[1]*w[0]), sp.expand(v[0]*w[2]-v[2]*w[0]), sp.expand(v[1]*w[2]-v[2]*w[1])]
    sol=sp.linsolve(eqs,b)
    for sl in sol:
        if any(x!=0 for x in sl):
            found+=1; print("polarity found:", sp.simplify(Bm.subs(dict(zip(b,sl)))).tolist()); break
    if found: break
print("side-swapping automorphisms tried:",tried," polarity exists:",found>0)
