"""Reproduce the frozen IPEDS study using Python 3.12, standard library only."""
import csv
import hashlib
import io
import json
import math
import platform
import random
import zipfile
from collections import defaultdict
from pathlib import Path
from statistics import median
from audit import Group, dot, envelope, independent

DATA = Path("data")
OUT = Path("results")
FILES = {
    2022: ("31cbfa5e68289525bff60e65b1094c610473ecc02d601461ae3875f854e2fde7", "adm2022_rv.csv"),
    2023: ("670ecc7c4313f044dcdf740a42ef10c6d2927eb9f890fb405e3569c8442aca27", "adm2023.csv"),
}
FIELDS = ("APPLCNW", "APPLCNM", "ADMSSNW", "ADMSSNM")


def panels(year):
    path = DATA / f"ADM{year}.zip"
    checksum, member = FILES[year]
    if hashlib.sha256(path.read_bytes()).hexdigest() != checksum:
        raise ValueError(f"Unexpected input: {path}")
    bins = defaultdict(list)
    flow = defaultdict(int)
    with zipfile.ZipFile(path) as z, z.open(member) as raw:
        for row in csv.DictReader(io.TextIOWrapper(raw, encoding="utf-8-sig")):
            flow["raw"] += 1
            if any(row["X"+f].strip() != "R" for f in FIELDS):
                flow["nonreported"] += 1
                continue
            numbers = tuple(int(row[f]) for f in FIELDS)
            w, m, yw, ym = numbers
            if min(w,m) <= 0 or not (0 <= yw <= w and 0 <= ym <= m):
                flow["invalid"] += 1
                continue
            if w+m < 100:
                flow["below_100"] += 1
                continue
            size = "small" if w+m < 1000 else "medium" if w+m < 10000 else "large"
            bins[size].append((row["UNITID"], *numbers))
            flow[size] += 1
    selected = []
    for size, rows in sorted(bins.items()):
        rows.sort(key=lambda row: hashlib.sha256(("P7:"+row[0]).encode()).digest())
        for start in range(0, min(len(rows)//8*8,80), 8):
            selected.append((size, start//8, rows[start:start+8]))
    return selected, dict(flow)


def rounded(y, n, step):
    # Integer inversion, including the half-up tie convention.
    middle = step*((2*y+step)//(2*step))
    return max(0,middle-step//2), min(n,middle+(step-1)//2)


def privatized(y, n, rng, radius):
    # Public test data only. Each outcome flip changes one count by one.
    u = rng.random()-0.5
    noise = -math.copysign(1, u)*math.log1p(-2*abs(u))
    z = min(n, max(0, y+noise))
    low, high = max(0,math.ceil(z-radius)), min(n,math.floor(z+radius))
    return low, high, z


def make_case(rows, step=1, rng=None):
    totals = [w+m for _,w,m,_,_ in rows]
    anchor = tuple(n/sum(totals) for n in totals)
    groups, plugin = [], []
    for g in (0,1):
        lo, hi, q, point = [], [], [], []
        for row, total in zip(rows,totals):
            n, y = row[1+g], row[3+g]
            if rng is None:
                l, h = rounded(y,n,step)
                z = (l+h)/2
            else:
                l,h,z = privatized(y,n,rng,math.log(2*len(rows)/0.05))
            q.append(n/total)
            lo.append(l/n)
            hi.append(h/n)
            point.append(z/n)
        groups.append(Group(tuple(q),tuple(lo),tuple(hi)))
        plugin.append(Group(tuple(q),tuple(point),tuple(point)))
    return anchor,groups,plugin


def centered(anchor, groups, rho):
    # Midpoint plug-in plus a uniform bound from conditional TV distance.
    centers = []
    radius = 0
    for g in groups:
        r = tuple((h+l)/2 for h,l in zip(g.high,g.low))
        centers.append(Group(g.share,r,r))
        peak = max(g.share)
        mix = rho*peak/((1-rho)*dot(anchor,g.share)+rho*peak)
        radius += mix*max(h-l for h,l in zip(g.high,g.low))
    e = envelope(anchor,centers,rho)
    return e.lower.value-radius, e.upper.value+radius


def write_csv(name, rows):
    with (OUT/name).open("w",newline="",encoding="utf-8") as f:
        writer = csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run():
    OUT.mkdir(exist_ok=True)
    profiles, simulations, units, flows = [], [], [], {}
    for year in FILES:
        selected, flows[year] = panels(year)
        for size, panel, rows in selected:
            key = dict(year=year,size=size,panel=panel)
            units.extend(dict(**key,unitid=r[0],applications_w=r[1],applications_m=r[2],
                              admitted_w=r[3],admitted_m=r[4]) for r in rows)
            for step in (1,10,100):
                anchor, groups, _ = make_case(rows,step)
                for rho in (0,0.01,0.10,0.25,1):
                    e = envelope(anchor,groups,rho)
                    il,iu = independent(anchor,groups,rho)
                    cl,cu = centered(anchor,groups,rho)
                    profiles.append(dict(**key,step=step,rho=rho,lo=e.lower.value,hi=e.upper.value,
                                single_lo=e.single_lower,single_hi=e.single_upper,
                                independent_lo=il,independent_hi=iu,centered_lo=cl,centered_hi=cu,
                                low_i=rows[e.lower.first][0],low_j=rows[e.lower.second][0],
                                low_t=e.lower.fraction_first,high_i=rows[e.upper.first][0],
                                high_j=rows[e.upper.second][0],high_t=e.upper.fraction_first))
            if panel == 0:
                anchor, exact, _ = make_case(rows)
                truth = envelope(anchor,exact,0.1)
                for rep in range(100):
                    rng = random.Random(f"P7-Laplace:{year}:{size}:{rep}")
                    anchor,box,plugin = make_case(rows,rng=rng)
                    e, p = envelope(anchor,box,.1), envelope(anchor,plugin,.1)
                    il,iu = independent(anchor,box,.1)
                    simulations.append(dict(**key,rep=rep,truth_lo=truth.lower.value,
                        truth_hi=truth.upper.value,lo=e.lower.value,hi=e.upper.value,
                        plugin_lo=p.lower.value,plugin_hi=p.upper.value,independent_lo=il,
                        independent_hi=iu,box_covers=int(all(l<=v<=h for g,x in zip(box,exact)
                        for l,v,h in zip(g.low,x.low,g.high)))))
    write_csv("institutions.csv",units)
    write_csv("profiles.csv",profiles)
    write_csv("private_releases.csv",simulations)
    summary = dict(platform=platform.platform(),python=platform.python_version(),flow=flows,
                   panels=len(units)//8,unique_institutions=len({r["unitid"] for r in units}),
                   primary=[],private=[])
    for size in ("small","medium","large"):
        for step in (1,10,100):
            sample = [r for r in profiles if r["size"]==size and r["step"]==step and r["rho"]==.1]
            def width(r,method):
                prefix = method+"_" if method else ""
                return 100*(r[prefix+"hi"]-r[prefix+"lo"])
            widths = {k or "sharp":median(width(r,k) for r in sample)
                      for k in ("","independent","centered","single")}
            summary["primary"].append(dict(size=size,step=step,n=len(sample),width_pp=widths,
                interior_count=sum(max(r["single_lo"]-r["lo"],r["hi"]-r["single_hi"])>1e-9 for r in sample)))
        sim = [r for r in simulations if r["size"]==size]
        def coverage(pre):
            return sum(r[pre+"lo"]<=r["truth_lo"]+1e-12 and r[pre+"hi"]>=r["truth_hi"]-1e-12 for r in sim)
        summary["private"].append(dict(size=size,n=len(sim),box=sum(r["box_covers"] for r in sim),
              shared=coverage(""),plugin=coverage("plugin_"),independent=coverage("independent_"),
              median_sharp_width_pp=median(100*(r["hi"]-r["lo"]) for r in sim),
              median_independent_width_pp=median(100*(r["independent_hi"]-r["independent_lo"]) for r in sim)))
    (OUT/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    run()
