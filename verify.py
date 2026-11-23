"""Verify returned witnesses, method nesting, and a deterministic full replay."""
import csv
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path
import replay
import study

ROOT = Path(__file__).resolve().parent
OUT = ROOT/"results"
NAMES = ("institutions.csv","profiles.csv","private_releases.csv","summary.json","next_year.csv","replay.json")


def read(name):
    with (OUT/name).open(encoding="utf-8") as f:
        return list(csv.DictReader(f))


def attained(anchor, groups, u, rho, orientation):
    # Direct conditional site probabilities, independent of pair solver.
    w = [(1-rho)*b+rho*t for b,t in zip(anchor,u)]
    result = 0
    for g,s in zip(groups,(1,-1)):
        d0 = sum(b*q for b,q in zip(anchor,g.share))
        d1 = sum(v*q for v,q in zip(w,g.share))
        c = [s*q*(v/d1-b/d0) for v,b,q in zip(w,anchor,g.share)]
        result += sum(a*(h if orientation*a>0 else l) for a,l,h in zip(c,g.low,g.high))
    return result


def main():
    subprocess.run([sys.executable,"-m","unittest","discover","-s","tests","-v"],check=True)
    cases = {(year,size,idx):rows for year in study.FILES for size,idx,rows in study.panels(year)[0]}
    profiles = read("profiles.csv")
    exact = {}
    worst_error = 0
    for r in profiles:
        key = int(r["year"]),r["size"],int(r["panel"])
        rows = cases[key]
        step,rho = int(r["step"]),float(r["rho"])
        anchor,groups,_ = study.make_case(rows,step)
        ids = [z[0] for z in rows]
        lo,hi = float(r["lo"]),float(r["hi"])
        for label,expected,orientation in (("low",lo,-1),("high",hi,1)):
            u = [0]*len(rows)
            t = float(r[label+"_t"])
            if not 0<=t<=1:
                raise ValueError("invalid mixture")
            u[ids.index(r[label+"_i"])] += t
            u[ids.index(r[label+"_j"])] += 1-t
            worst_error = max(worst_error,abs(attained(anchor,groups,u,rho,orientation)-expected))
        if step == 1:
            exact[key,rho] = (lo,hi)
        truth_lo,truth_hi = exact[key,rho]
        if not (lo<=truth_lo+1e-12 and truth_hi<=hi+1e-12 and lo<=1e-12 and hi>=-1e-12
                and float(r["single_lo"])>=lo-1e-12 and float(r["single_hi"])<=hi+1e-12
                and all(float(r[p+"lo"])<=lo+1e-12 and float(r[p+"hi"])>=hi-1e-12
                        for p in ("independent_","centered_"))):
            raise AssertionError("profile containment")
    if len(profiles)!=900 or len(read("private_releases.csv"))!=600 or worst_error>1e-12:
        raise AssertionError("run size or witness verification failed")
    # Repeat in isolation. CSVs and replay summary must match byte for byte.
    # summary.json includes the platform/Python version, compared semantically
    # so an independent user may also reproduce on another OS.
    from contextlib import redirect_stdout
    import io
    with tempfile.TemporaryDirectory() as temp:
        study.OUT = replay.OUT = Path(temp)
        with redirect_stdout(io.StringIO()):
            study.run()
            replay.run()
        for name in NAMES:
            if name.endswith(".json"):
                expected,actual = [json.loads((p/name).read_text()) for p in (OUT,Path(temp))]
                if name=="summary.json":
                    for obj in (expected,actual):
                        for k in ("platform","python"):
                            obj.pop(k)
                assert expected==actual,name
            else:
                assert (OUT/name).read_bytes()==(Path(temp)/name).read_bytes(),name
    print(f"Verified {len(profiles)*2} witnesses; max absolute error {worst_error:.3g}.")
    print("All six result files reproduced.")
    for name in NAMES:
        print(hashlib.sha256((OUT/name).read_bytes()).hexdigest(),name)


if __name__=="__main__":
    main()
