"""Secondary performance of shared-table screening and directional decisions."""
import csv
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from statistics import median
from audit import independent, linear_certificate
from edfacts_audit import load, projections

OUT = Path("results")
FILE = "edfacts_screen.json"


def rational_fixed_case(sites):
    """Exact rational support at the vertices, independently of the integer sweep."""
    total = sum(s.total for s in sites)
    base = [Fraction(s.total,total) for s in sites]
    target = [Fraction(1,len(sites)) for s in sites]
    coeff = []
    for counts in ([s.women for s in sites],[s.men for s in sites]):
        base_mass = sum(b*n/s.total for b,n,s in zip(base,counts,sites))
        target_mass = sum(w*n/s.total for w,n,s in zip(target,counts,sites))
        coeff.append([(w/target_mass-b/base_mass)/s.total for w,b,s in zip(target,base,sites)])
    limits = [Fraction(0) for _ in range(4)]
    for s,c,d in zip(sites,*coeff):
        v = s.vertices()
        paired = [c*x-d*y for x,y in v]
        xs,ys = zip(*v)
        marginal = [c*x-d*y for x in (min(xs),max(xs)) for y in (min(ys),max(ys))]
        for i,value in enumerate((min(paired),max(paired),min(marginal),max(marginal))):
            limits[i] += value
    return limits


def analyze():
    with (OUT/"edfacts_districts.csv").open(encoding="utf-8") as f:
        rows=list(csv.DictReader(f))
    selected,_=load(False)
    primary=[r for r in rows if r["rule"]=="nearest" and r["setting"]=="0.10"]
    tally=dict(separate_pass=0,scan_pass=0,scan_violation=0,scan_undecided=0,
               exact_pass=0,undecided_exact_pass=0)
    widths=[[],[],[]]
    for row in primary:
        _,sites=selected[row["district"]]
        n=sum(s.total for s in sites)
        anchor=tuple(s.total/n for s in sites)
        groups=projections(sites)
        # Each extreme of a level gap uses opposite-coordinate corners.
        sl,su=independent(anchor,groups,.1)
        bracket=linear_certificate(anchor,groups,.1)
        lo,hi=float(row["joint_lo"]),float(row["joint_hi"])
        assert sl-1e-10<=lo<=hi<=su+1e-10
        assert bracket.lower-1e-10<=lo<=bracket.inner_lower+1e-10
        assert bracket.inner_upper-1e-10<=hi<=bracket.upper+1e-10
        passed=max(-lo,hi)<=.05+1e-12
        tally["exact_pass"]+=passed
        tally["separate_pass"]+=max(-sl,su)<=.05+1e-12
        if max(-bracket.lower,bracket.upper)<=.05+1e-12:
            assert passed
            tally["scan_pass"]+=1
        elif max(-bracket.inner_lower,bracket.inner_upper)>.05+1e-12:
            assert not passed
            tally["scan_violation"]+=1
        else:
            tally["scan_undecided"]+=1
            tally["undecided_exact_pass"]+=passed
        for values,width in zip(widths,(su-sl,bracket.upper-bracket.lower,hi-lo)):
            values.append(100*width)
    directions={}
    witness_checks={}
    def sign(lo,hi):
        return 1 if lo>1e-12 else -1 if hi< -1e-12 else 0
    for rule in ("nearest","union"):
        cases=[r for r in rows if r["rule"]==rule and r["setting"]=="equal_school"]
        changed=[]
        for r in cases:
            joint=sign(float(r["joint_lo"]),float(r["joint_hi"]))
            proj=sign(float(r["proj_lo"]),float(r["proj_hi"]))
            if joint!=proj:
                changed.append(dict(district=r["district"],name=r["name"],joint=joint,projected=proj,
                             joint_pp=[100*float(r[k]) for k in ("joint_lo","joint_hi")],
                             projection_pp=[100*float(r[k]) for k in ("proj_lo","proj_hi")]))
        directions[rule]=changed
        # Post hoc check of the one fragile sign difference from the nearest run.
        district = "0628320"
        case = next(r for r in cases if r["district"]==district)
        sites = selected[district][1] if rule=="nearest" else load(True)[0][district][1]
        exact = rational_fixed_case(sites)
        fields = ("joint_lo","joint_hi","proj_lo","proj_hi")
        for k,v in zip(fields,exact):
            assert abs(float(case[k])-float(v))<1e-12
        witness_checks[rule] = dict(zip(fields,map(str,exact)))
    return dict(primary_districts=len(primary),decisions=tally,
                median_width_pp=dict(zip(("separate","scan","joint"),map(median,widths))),direction_changes=directions,
                old_adobe_rational_change=witness_checks)


if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify",action="store_true")
    args=parser.parse_args()
    data=(json.dumps(analyze(),indent=2)+"\n").encode()
    if args.verify:
        assert data==(OUT/FILE).read_bytes()
    else:
        (OUT/FILE).write_bytes(data)
    print(data.decode(),hashlib.sha256(data).hexdigest(),FILE)
