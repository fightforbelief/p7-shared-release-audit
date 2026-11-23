"""Paired comparisons on previously frozen EDFacts inputs and decisions."""
import argparse
import csv
import hashlib
import io
import json
from collections import defaultdict
from itertools import groupby
from pathlib import Path
from statistics import median
from time import perf_counter

from audit import envelope, linear_certificate
from edfacts_audit import load, projections
from shared_baselines import METHODS, bounds
from totals_audit import endpoints

OUT = Path("results")
NAMES = ("chord",*METHODS)
EPS = 1e-12


def inputs(sites):
    total = sum(s.total for s in sites)
    return tuple(s.total/total for s in sites),projections(sites)


def radius(interval):
    return max(-interval[0],interval[1])


def decision(inner,outer,tol):
    if radius(inner)>tol+EPS:
        return "violation"
    if radius(outer)<=tol+EPS:
        return "pass"
    return "open"


def paired():
    truth = {}
    with (OUT/"edfacts_districts.csv").open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row["setting"]!="equal_school":
                truth[row["rule"],row["setting"],row["district"]] = row
    rows = []
    for rule,subset in groupby(sorted(truth),key=lambda key:key[0]):
        selected,_=load(rule=="union")
        for _,setting,district in subset:
            sites = selected[district][1]
            anchor,groups = inputs(sites)
            rho = float(setting)
            vals = bounds(anchor,groups,rho)
            br = linear_certificate(anchor,groups,rho)
            vals["chord"] = (br.lower,br.upper)
            exact = truth[rule,setting,district]
            jl,ju,pl,pu = (float(exact[key]) for key in ("joint_lo","joint_hi","proj_lo","proj_hi"))
            il,iu = vals["inner"]
            assert abs(il-br.inner_lower)<EPS and abs(iu-br.inner_upper)<EPS
            assert jl<=il+EPS and iu<=ju+EPS
            for name in NAMES:
                lo,hi = vals[name]
                assert lo<=pl+EPS and pu<=hi+EPS
                row = dict(rule=rule,setting=setting,district=district,method=name,
                           inner_lo=il,inner_hi=iu,outer_lo=lo,outer_hi=hi)
                for pp in (2.5,5,7.5):
                    d = decision((il,iu),(lo,hi),pp/100)
                    passes = radius((jl,ju))<=pp/100+EPS
                    assert d=="open" or (d=="pass")==passes
                    row[str(pp)]=d
                rows.append(row)
    summary = []
    by_case = defaultdict(dict)
    for r in rows:
        by_case[r["rule"],r["setting"],r["district"]][r["method"]]=r
    for (rule,setting),cases in groupby(sorted(by_case),key=lambda key:key[:2]):
        districts = [by_case[key] for key in cases]
        for pp in ("2.5","5","7.5"):
            for method in NAMES:
                decisions = [d[method][pp] for d in districts]
                scans = [d["chord"][pp] for d in districts]
                summary.append(dict(rule=rule,setting=setting,tolerance_pp=float(pp),method=method,
                    passed=decisions.count("pass"),violation=decisions.count("violation"),
                    unresolved=decisions.count("open"),
                    chord_only=sum(s!="open" and d=="open" for s,d in zip(scans,decisions)),
                    comparator_only=sum(d!="open" and s=="open" for s,d in zip(scans,decisions)),
                    median_outer_width_pp=median(100*(d[method]["outer_hi"]-d[method]["outer_lo"]) for d in districts)))
    return rows,summary


def output(verify):
    for filename,rows in zip(("baseline_districts.csv","baseline_comparison.csv"),paired()):
        buf=io.StringIO(newline="")
        writer=csv.DictWriter(buf,fieldnames=list(rows[0]),lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        data=buf.getvalue().encode()
        path=OUT/filename
        if verify:
            assert path.read_bytes()==data,filename
        else:
            path.write_bytes(data)
        print(len(rows),hashlib.sha256(data).hexdigest(),filename)
        if filename=="baseline_comparison.csv":
            for r in rows:
                if (r["rule"],r["setting"],r["tolerance_pp"])==("nearest","0.10",5):
                    print(r)


def benchmark():
    selected,_=load(False)
    cases=[(*inputs(sites),sites) for _,(_,sites) in sorted(selected.items())]

    def run(method, fallback):
        start=perf_counter()
        counts=[0,0,0]  # within tolerance, witnessed violation, projected fallbacks
        for anchor,groups,sites in cases:
            d="open"
            if method=="chord":
                br=linear_certificate(anchor,groups,.1)
                d=decision((br.inner_lower,br.inner_upper),(br.lower,br.upper),.05)
            else:
                v=bounds(anchor,groups,.1,(method,))
                d=decision(v["inner"],v[method],.05)
            if d=="open" and fallback:
                counts[2]+=1
                e=envelope(anchor,groups,.1)
                if radius((e.lower.value,e.upper.value))<=.05+EPS:
                    d="pass"
                else:
                    lo,hi=endpoints(anchor,sites,.1)
                    d="pass" if radius((lo.value,hi.value))<=.05+EPS else "violation"
            if d=="pass":
                counts[0]+=1
            elif d=="violation":
                counts[1]+=1
        return perf_counter()-start,counts

    # Direct primary projected and joint radii were previously verified equal.
    # Time direct projection alone as the strongest exact reference on this cohort.
    def direct():
        start=perf_counter()
        passed=0
        for anchor,groups,_ in cases:
            e=envelope(anchor,groups,.1)
            passed+=radius((e.lower.value,e.upper.value))<=.05+EPS
        return perf_counter()-start,[passed,len(cases)-passed,len(cases)]

    methods=(*NAMES,"direct")
    records={m:[] for m in methods}
    for rep in range(6):
        for m in methods[rep:]+methods[:rep]:
            sec,counts=direct() if m=="direct" else run(m,True)
            scan_sec,scan_counts=(None,None) if m=="direct" else run(m,False)
            assert counts[:2]==[192,225]
            if rep:
                records[m].append(dict(seconds=sec,screen_only_seconds=scan_sec,
                                       first_pass_counts=scan_counts,projected_calls=counts[2]))
    report=dict(scope="417 prepared inputs; standalone screen followed by exact fallback; direct projected reference",
                warm_runs=5,methods=records)
    (OUT/"baseline_timing.json").write_bytes((json.dumps(report,indent=2)+"\n").encode())
    for name,times in records.items():
        print(name,median(r["seconds"] for r in times),times[0]["projected_calls"])


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--verify",action="store_true")
    mode.add_argument("--timing",action="store_true")
    args=parser.parse_args()
    if args.timing:
        benchmark()
    else:
        output(args.verify)
