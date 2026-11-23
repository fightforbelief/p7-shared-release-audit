"""Release construction, paired audit, and repeatability checks for public totals."""
import argparse
import csv
import hashlib
import io
import json
import tempfile
import zipfile
from pathlib import Path
from statistics import median
from audit import envelope
from study import DATA, FILES, panels, rounded
from totals_audit import Site, endpoints, linear, projections

OUT = Path("results")
RHO, LIMIT = .1,.05
CASES = ("all_groups_rounded","residual_unreleased")
FILES_OUT = ("public_total_panels.csv","public_total_summary.json")


def read_csv(path):
    with path.open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))


def certify(lo,hi):
    return max(-lo,hi) <= LIMIT


def number(row,field):
    value = row[field].strip()
    return int(value) if value else None


def site_from_release(row, use_residual):
    nw,nm = (int(row["APPLCN"+s]) for s in ("W","M"))
    capacity = int(row["APPLCN"])-nw-nm
    admitted = int(row["ADMSSN"])
    if capacity < 0:
        raise ValueError("negative residual application capacity")
    low,high = 0,capacity
    if use_residual:
        caps = [number(row,"APPLCN"+s) for s in ("AN","UN")]
        if sum(c for c in caps if c is not None)>capacity:
            raise ValueError("residual category denominators exceed capacity")
        lows,highs = [],[]
        for i,suffix in enumerate(("AN","UN")):
            cap = caps[i] if caps[i] is not None else capacity-(caps[1-i] or 0)
            value = number(row,"ADMSSN"+suffix)
            l,h = (0,cap) if value is None else rounded(value,cap,10)
            lows.append(l)
            highs.append(h)
        low,high = sum(lows),min(capacity,sum(highs))
    sum_limits = max(0,admitted-high),min(nw+nm,admitted-low)
    return Site(nw,nm,rounded(int(row["ADMSSNW"]),nw,10),
                rounded(int(row["ADMSSNM"]),nm,10),sum_limits)


def direct(anchor,sites,rho,witness,direction):
    # Enumerate integer pairs in the original cells at the returned weights.
    u = [0.0]*len(sites)
    u[witness.first] += witness.fraction_first
    u[witness.second] += 1-witness.fraction_first
    target = [(1-rho)*b+rho*v for b,v in zip(anchor,u)]
    def masses(w):
        return [sum(v*n/s.total for v,s,n in zip(w,sites,ns))
                for ns in ([s.women for s in sites],[s.men for s in sites])]
    bw,bm = masses(anchor)
    dw,dm = masses(target)
    extremum = min if direction==-1 else max
    total = 0
    for b,w,s in zip(anchor,target,sites):
        cw,cm = (w/dw-b/bw)/s.total,(b/bm-w/dm)/s.total
        total += extremum(cw*x+cm*y for x in range(s.women_bounds[0],s.women_bounds[1]+1)
                          for y in range(s.men_bounds[0],s.men_bounds[1]+1) if s.contains(x,y))
    return total


def build(out):
    original = {(int(r["year"]),r["size"],int(r["panel"])):r for r in read_csv(OUT/"profiles.csv")
                if float(r["rho"])==RHO and int(r["step"])==10}
    old_linear = {(int(r["year"]),r["size"],int(r["panel"])):r for r in read_csv(OUT/"curvature.csv")
                  if float(r["rho"])==RHO and int(r["step"])==10}
    records = []
    max_error = 0
    for year in FILES:
        selected,_ = panels(year)  # Verifies official archive checksum.
        with zipfile.ZipFile(DATA/f"ADM{year}.zip") as z, z.open(FILES[year][1]) as raw:
            full = {r["UNITID"]:r for r in csv.DictReader(io.TextIOWrapper(raw,encoding="utf-8-sig"))}
        for size,idx,panel in selected:
            key = year,size,idx
            rows = [full[r[0]] for r in panel]
            assert all(r["XADMSSN"].strip()==r["XAPPLCN"].strip()=="R" for r in rows)
            n = [r[1]+r[2] for r in panel]
            b = tuple(v/sum(n) for v in n)
            old = original[key]
            bl,bu = float(old["lo"]),float(old["hi"])
            old_scan = old_linear[key]
            box_linear = certify(float(old_scan["lo"]),float(old_scan["hi"]))
            for name in CASES:
                sites = [site_from_release(r,name==CASES[0]) for r in rows]
                # Ground truth is used only for validation, never to tighten the band.
                assert all(s.contains(p[3],p[4]) for s,p in zip(sites,panel))
                lo,hi = endpoints(b,sites,RHO)
                scan = linear(b,sites,RHO)
                rectangular = envelope(b,projections(sites),RHO)
                eps = 1e-12
                assert bl-eps<=rectangular.lower.value<=lo.value+eps
                assert hi.value-eps<=rectangular.upper.value<=bu+eps
                assert scan.lower-eps<=lo.value<=scan.inner_lower+eps
                assert scan.inner_upper-eps<=hi.value<=scan.upper+eps
                for witness,sign in ((lo,-1),(hi,1)):
                    max_error = max(max_error,abs(direct(b,sites,RHO,witness,sign)-witness.value))
                records.append(dict(year=year,size=size,panel=idx,release=name,
                    box_lo=bl,box_hi=bu,box_exact=certify(bl,bu),box_linear=box_linear,
                    lo=lo.value,hi=hi.value,scan_lo=scan.lower,scan_hi=scan.upper,
                    exact=certify(lo.value,hi.value),linear=certify(scan.lower,scan.upper),
                    projected_lo=rectangular.lower.value,projected_hi=rectangular.upper.value,
                    low_i=panel[lo.first][0],low_j=panel[lo.second][0],low_t=lo.fraction_first,
                    high_i=panel[hi.first][0],high_j=panel[hi.second][0],high_t=hi.fraction_first))
    assert len(records)==120 and max_error<1e-12
    for a,b in zip(records[::2],records[1::2]):
        assert b["lo"]-1e-12<=a["lo"]<=a["hi"]<=b["hi"]+1e-12
    summary = []
    for name in CASES:
        for size in ("small","medium","large","all"):
            rows = [r for r in records if r["release"]==name and size in ("all",r["size"])]
            yes = [r for r in rows if r["box_exact"]]
            def count(rs,k): return sum(r[k] for r in rs)
            summary.append(dict(release=name,size=size,n=len(rows),old_exact=len(yes),
                old_linear=count(rows,"box_linear"),retained=count(yes,"exact"),
                linear_within_old=count(yes,"linear"),exact=count(rows,"exact"),linear=count(rows,"linear"),
                added=[f"{r['year']}/{r['panel']}/{r['size']}" for r in rows if r["exact"] and not r["box_exact"]],
                exact_only=[f"{r['year']}/{r['panel']}/{r['size']}" for r in rows if r["exact"] and not r["linear"]],
                median_width_pp=median(100*(r["hi"]-r["lo"]) for r in rows),
                paired_median_reduction_pp=median(100*(r["box_hi"]-r["box_lo"]-r["hi"]+r["lo"]) for r in rows),
                max_projection_endpoint_error_pp=max(100*max(r["lo"]-r["projected_lo"],r["projected_hi"]-r["hi"]) for r in rows)))
    assert len([r for r in records if r["release"]==CASES[0] and r["size"]=="small" and r["box_exact"]])==13
    with (out/FILES_OUT[0]).open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]))
        w.writeheader()
        w.writerows(records)
    (out/FILES_OUT[1]).write_bytes((json.dumps(dict(max_witness_error=max_error,groups=summary),indent=2)+"\n").encode())
    return summary


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify",action="store_true")
    args = parser.parse_args()
    if args.verify:
        with tempfile.TemporaryDirectory() as d:
            build(Path(d))
            for file in FILES_OUT:
                data = (OUT/file).read_bytes()
                assert data==(Path(d)/file).read_bytes(),file
                print(hashlib.sha256(data).hexdigest(),file)
    else:
        for item in build(OUT):
            print(json.dumps(item))
