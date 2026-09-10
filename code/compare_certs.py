import json, sys, itertools
sys.path.insert(0,'/home/claude/repo/code')
from classify_V_general import Layout, V, E, S1, S2, S3, certificate, vmul
from solve_general import dotpoly, poly_str

def old_coords(lay):
    pts={0:((1,None),(0,None),(0,None))}
    pts[lay.T(1,0)]=((0,None),(1,None),(1,'t0')); pts[lay.T(1,1)]=((0,None),(1,None),(-1,'t0'))
    pts[lay.T(2,0)]=((0,None),(1,None),(1,'t1')); pts[lay.T(2,1)]=((0,None),(1,None),(-1,'t1'))
    if lay.b_axis==3:
        pts[lay.B(0)]=((1,None),(1,'S0'),(0,None)); pts[lay.B(1)]=((1,None),(-1,'S0'),(0,None))
    else:
        pts[lay.B(0)]=((1,None),(0,None),(1,'S0')); pts[lay.B(1)]=((1,None),(0,None),(-1,'S0'))
    for j in range(1,lay.c+1):
        for al in V:
            u=None if j==1 else 'u%d'%j; v=None if j==1 else 'v%d'%j
            pts[lay.P(j,al)]=((1,None),(al[0]*al[1],u),(al[0]*al[2],v))
    lns={0:((1,None),(0,None),(0,None))}
    lns[lay.Tp(1,0)]=((0,None),(1,None),(-1,None)); lns[lay.Tp(1,1)]=((0,None),(1,None),(1,None))
    lns[lay.Tp(2,0)]=((0,None),(1,'t1'),(-1,None)); lns[lay.Tp(2,1)]=((0,None),(1,'t1'),(1,None))
    for k,cen in enumerate(lay.bp):
        var='T%d'%k
        if cen==3:
            lns[lay.Bp(k,0)]=((1,None),(1,var),(0,None)); lns[lay.Bp(k,1)]=((1,None),(-1,var),(0,None))
        else:
            lns[lay.Bp(k,0)]=((1,None),(0,None),(1,var)); lns[lay.Bp(k,1)]=((1,None),(0,None),(-1,var))
    for a in range(1,lay.cp+1):
        for al in V:
            lns[lay.G(a,al)]=((1,None),(al[0]*al[1],'m%d'%a),(al[0]*al[2],'n%d'%a))
    return pts,lns

def decode(lay, gens):
    pts,lns=old_coords(lay)
    table={}
    for p in range(23):
        for l in range(23):
            poly=dotpoly(pts[p],lns[l])
            if poly is None: continue
            table.setdefault(poly_str(poly),[]).append((p,l))
    pl=[0]*23
    for p in range(23):
        for l in range(23):
            if dotpoly(pts[p],lns[l]) is None: pl[p]|=1<<l
    bad=[]
    for g in gens:
        g=g.replace(' ','')
        if g not in table: bad.append(g); continue
        for p,l in table[g]: pl[p]|=1<<l
    return pl,bad

def check(pl):
    lp=[0]*23
    for p in range(23):
        m=pl[p]
        while m:
            l=(m&-m).bit_length()-1; lp[l]|=1<<p; m&=m-1
    deg_ok=all(bin(x).count('1')==4 for x in pl) and all(bin(x).count('1')==4 for x in lp)
    c4=all(bin(pl[p]&pl[q]).count('1')<=1 for p,q in itertools.combinations(range(23),2))
    return deg_ok,c4

for tag,fname,cp,bax,bp in [('3_3_233','23_B_4_3_3_233',3,3,[2,3,3]),('3_3_223','23_B_4_3_3_223',3,3,[2,2,3])]:
    lay=Layout(4,cp,bax,bp)
    mine=set(json.load(open('/home/claude/all_%s.json'%tag)).keys())
    d=json.load(open('/home/claude/repo/certificates/%s_certificates.json'%fname))
    certs={}; nbad=0; nfail=0; keys_old=[]
    for idx,s in d['structures'].items():
        pl,bad=decode(lay,s['generators'])
        if bad: nbad+=1; continue
        ok1,ok2=check(pl)
        if not(ok1 and ok2): nfail+=1; continue
        ce=certificate(lay,pl).hex(); keys_old.append(ce)
        certs.setdefault(ce,[]).append(idx)
    olds=set(certs)
    dupes=sum(len(v)-1 for v in certs.values())
    print(f"{tag}: certified {len(d['structures'])}, undecodable {nbad}, failing degree/4-cycle {nfail}, "
          f"distinct up to isomorphism {len(olds)} (isomorphic duplicates {dupes})")
    print(f"      current enumeration {len(mine)};  current \\ certified: {len(mine-olds)};  certified \\ current: {len(olds-mine)}")
