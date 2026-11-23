"""Standard rational Bernstein subdivision on the unresolved California cases."""
import argparse
import csv
import hashlib
import json
from itertools import combinations
from pathlib import Path

from audit import lerp
from edfacts_audit import fixed_bounds, load, mixture
from edfacts_baselines import EPS, OUT, decision, inputs
from edfacts_extract import EXTRACT
from shared_baselines import VertexGroup
from totals_audit import endpoints


def upper(terms, pieces):
    best=0.
    for i,j in combinations(range(len(terms[0].den)),2):
        def coefficient(l,h):
            num=den=0.
            for g,other in ((terms[0],terms[1]),(terms[1],terms[0])):
                dl=lerp(other.den[j],other.den[i],l)
                dh=lerp(other.den[j],other.den[i],h)
                p=g.mass[i]*(l*dh+h*dl)
                q=g.mass[j]*((1-l)*dh+(1-h)*dl)
                num+=(p*(g.low[i]-g.base)+q*(g.low[j]-g.base)
                      +max(0.,p-(p+q)*g.share[i])*g.width[i]
                      +max(0.,q-(p+q)*g.share[j])*g.width[j])
                den+=lerp(g.den[j],g.den[i],l)*dh
            return num/den
        best=max(best,coefficient(0,0))
        for k in range(pieces):
            l,h=k/pieces,(k+1)/pieces
            best=max(best,coefficient(l,h),coefficient(h,h))
    return best


def results():
    with (OUT/"baseline_districts.csv").open() as f:
        remaining=[r for r in csv.DictReader(f) if r["method"]=="bernstein" and
                   any(r[key]=="open" for key in ("2.5","5","7.5"))]
    report=[]
    for r in remaining:
        district=r["district"]
        name,sites=load(r["rule"]=="union")[0][district]
        b,groups=inputs(sites)
        rho=float(r["setting"])
        with EXTRACT.open(encoding="utf-8") as f:
            ids=[row["NCESSCH"] for row in csv.DictReader(f)
                 if row["LEAID"]==district and row["CATEGORY"]=="ALL"]
        assert len(ids)==len(sites) and ids==sorted(ids)
        inn=tuple(float(r[key]) for key in ("inner_lo","inner_hi"))
        exact=endpoints(b,sites,rho)
        cases=[]
        for side,w in enumerate(exact):
            check=fixed_bounds(sites,b,mixture(b,rho,w))[side]
            assert abs(check-w.value)<EPS
            cases.append(dict(change_pp=100*w.value,school_i=ids[w.first],school_j=ids[w.second],
                              fraction_at_i=w.fraction_first))
        prev=(float(r["outer_lo"]),float(r["outer_hi"]))
        levels=[]
        for pieces in (1,2,4,8):
            vals=[]
            for direction in (-1,1):
                terms=[VertexGroup(b,g,s,rho) for g,s in zip(groups,(direction,-direction))]
                vals.append(direction*upper(terms,pieces))
            assert prev[0]-EPS<=vals[0]<=exact[0].value+EPS
            assert exact[1].value-EPS<=vals[1]<=prev[1]+EPS
            if pieces==1:
                assert all(abs(x-y)<EPS for x,y in zip(vals,prev))
            levels.append(dict(pieces=pieces,outer_pp=[100*x for x in vals],
                               decisions={v:decision(inn,vals,float(v)/100) for v in ("2.5","5","7.5")}))
            prev=vals
        report.append(dict(district=district,name=name,rho=rho,rule=r["rule"],exact=cases,levels=levels))
    return report


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify",action="store_true")
    args=parser.parse_args()
    data=(json.dumps(results(),indent=2)+"\n").encode()
    path=Path("results/bernstein_refinement.json")
    if args.verify:
        assert path.read_bytes()==data
    else:
        path.write_bytes(data)
    print(data.decode(),hashlib.sha256(data).hexdigest())
