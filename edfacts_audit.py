"""External audit of released ALL and LEP math proficiency bands."""
import argparse
import csv
import hashlib
import json
import re
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median
from audit import envelope
from edfacts_extract import EXTRACT
from totals_audit import Site, endpoints, projections

OUT = Path("results")
OUTPUTS = ("edfacts_districts.csv","edfacts_summary.json")
SOURCE_PIN = "42e47099976d5f3a8b85aa105a8423817e258ac94a3d45ac9f0abfe3943f323a"


def count_band(n, label, conservative=False):
    if label == "PS":
        return 0,n
    if label.isdigit():
        low = high = int(label)
    elif re.fullmatch(r"[0-9]+-[0-9]+",label):
        low,high = map(int,label.split("-"))
    else:
        match = re.fullmatch(r"(LT|LE|GT|GE)([0-9]+)",label)
        if not match:
            raise ValueError(f"unknown percent label {label!r}")
        op, v = match[1],int(match[2])
        low,high = {"LT":(0,v-1),"LE":(0,v),"GT":(v+1,100),"GE":(v,100)}[op]
    if not 0<=low<=high<=100:
        raise ValueError("invalid percent range")
    # Include both possible nearest-percent ties. The sensitivity also
    # admits a percentage truncated to an integer (strict upper edge).
    a = max(0,-(-n*(2*low-1)//200))
    b = min(n,n*(2*high+1)//200)
    if conservative:
        b = max(b,min(n,(n*(high+1)-1)//100))
    return a,b


def load(conservative):
    if hashlib.sha256(EXTRACT.read_bytes()).hexdigest()!=SOURCE_PIN:
        raise ValueError("public extract checksum mismatch")
    districts = defaultdict(lambda:defaultdict(dict))
    names = {}
    with EXTRACT.open(encoding="utf-8",newline="") as f:
        for r in csv.DictReader(f):
            district = r["LEAID"]
            row = districts[district][r["NCESSCH"]]
            if r["CATEGORY"] in row:
                raise ValueError("duplicate school-category")
            row[r["CATEGORY"]] = r
            names[district] = r["LEANM"]
    flow = Counter()
    selected = {}
    for district,schools in districts.items():
        sites = []
        incomplete = False
        had_all = False
        for rows in schools.values():
            if "ALL" not in rows:
                flow["orphan_subgroup_rows"] += 1
                continue
            had_all = True
            flow["schools_all"] += 1
            if any(k not in rows or not rows[k]["NUMVALID"].isdigit() or
                   rows[k]["PCTPROF"] in ("",".") for k in ("LEP","ALL")):
                flow["missing_row_count_or_percent"] += 1
                incomplete = True
                continue
            eg,allrow = rows["LEP"],rows["ALL"]
            e,n = int(eg["NUMVALID"]),int(allrow["NUMVALID"])
            if not 0 < e < n:
                flow["zero_or_invalid_group"] += 1
                incomplete = True
                continue
            l,h = count_band(e,eg["PCTPROF"],conservative)
            a,b = count_band(n,allrow["PCTPROF"],conservative)
            if l>h or a>b or h+n-e<a or l>b:
                flow["empty_polygon"] += 1
                incomplete = True
                continue
            sites.append(Site(e,n-e,(l,h),(0,n-e),(a,b)))
        if not had_all:
            continue
        flow["districts_all"] += 1
        if incomplete:
            flow["incomplete_districts"] += 1
        elif not 3<=len(sites)<=25:
            flow["size_excluded_districts"] += 1
        else:
            selected[district] = (names[district],sites)
            flow["selected_schools"] += len(sites)
    flow["selected_districts"] = len(selected)
    return selected,dict(sorted(flow.items()))


def fixed_bounds(sites,b,w):
    # Independent check via a one-dimensional integer sweep. For each
    # feasible learner count, the complement extrema occur at its two limits.
    def mass(weights):
        return (sum(z*s.women/s.total for z,s in zip(weights,sites)),
                sum(z*s.men/s.total for z,s in zip(weights,sites)))
    b1,b0 = mass(b)
    w1,w0 = mass(w)
    jl=ju=pl=pu=0.
    for s,bs,ws in zip(sites,b,w):
        c = (ws/w1-bs/b1)/s.total
        d = (bs/b0-ws/w0)/s.total
        vertices = s.vertices()
        xx,yy = zip(*vertices)
        xl,xh,yl,yh = min(xx),max(xx),min(yy),max(yy)
        pl += min(c*xl,c*xh)+min(d*yl,d*yh)
        pu += max(c*xl,c*xh)+max(d*yl,d*yh)
        vals = []
        for x in range(xl,xh+1):
            low = max(yl,s.sum_bounds[0]-x)
            high = min(yh,s.sum_bounds[1]-x)
            vals.extend((c*x+d*low,c*x+d*high))
        jl += min(vals)
        ju += max(vals)
    return jl,ju,pl,pu


def mixture(b,rho,witness):
    w = [(1-rho)*v for v in b]
    w[witness.first] += rho*witness.fraction_first
    w[witness.second] += rho*(1-witness.fraction_first)
    return w


def audit(sites,setting):
    tested = sum(s.total for s in sites)
    b = [s.total/tested for s in sites]
    if setting=="equal_school":
        return fixed_bounds(sites,b,[1/len(sites)]*len(sites))
    rho = float(setting)
    lo,hi = endpoints(b,sites,rho)
    rect = envelope(b,projections(sites),rho)
    for idx,v in enumerate((lo,hi,rect.lower,rect.upper)):
        calc = fixed_bounds(sites,b,mixture(b,rho,v))[idx]
        if abs(calc-v.value)>1e-11:
            raise AssertionError("witness failed direct verification")
    return lo.value,hi.value,rect.lower.value,rect.upper.value


def build(out):
    rows = []
    flows = {}
    for conservative in (False,True):
        rule = "union" if conservative else "nearest"
        selected,flows[rule] = load(conservative)
        for district,(name,sites) in sorted(selected.items()):
            settings = ("0.10","equal_school") if conservative else ("0.10","0.25","1.00","equal_school")
            for setting in settings:
                jl,ju,pl,pu = audit(sites,setting)
                if jl<pl-1e-11 or ju>pu+1e-11:
                    raise AssertionError("joint interval exceeds projection")
                j,p = max(-jl,ju),max(-pl,pu)
                rows.append(dict(rule=rule,district=district,name=name,schools=len(sites),
                                 setting=setting,joint_lo=jl,joint_hi=ju,proj_lo=pl,proj_hi=pu,
                                 joint_radius_pp=100*j,proj_radius_pp=100*p,
                                 radius_gain_pp=max(0.,100*(p-j)),endpoint_gain_pp=max(0.,100*(jl-pl),100*(pu-ju))))
    results = []
    for rule,setting in sorted({(r["rule"],r["setting"]) for r in rows}):
        group = [r for r in rows if (r["rule"],r["setting"])==(rule,setting)]
        largest = max(group,key=lambda r:r["radius_gain_pp"])
        decisions = {}
        for tau in (2.5,5,7.5):
            joints = sum(r["joint_radius_pp"]<=tau+1e-10 for r in group)
            projs = sum(r["proj_radius_pp"]<=tau+1e-10 for r in group)
            flips = [r for r in group if r["joint_radius_pp"]<=tau+1e-10<r["proj_radius_pp"]]
            decisions[str(tau)]=dict(joint=joints,projection=projs,flips=len(flips),
                material_flips=sum(r["radius_gain_pp"]>=.1 for r in flips))
        results.append(dict(rule=rule,setting=setting,n=len(group),decisions=decisions,
                unequal_radii=sum(r["radius_gain_pp"]>1e-8 for r in group),
                material_radius=sum(r["radius_gain_pp"]>=.1 for r in group),
                median_radius_gain_pp=median(r["radius_gain_pp"] for r in group),
                max_endpoint_gain_pp=max(r["endpoint_gain_pp"] for r in group),
                largest_radius_case=largest))
    with (out/OUTPUTS[0]).open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (out/OUTPUTS[1]).write_bytes((json.dumps(dict(flow=flows,comparisons=results),indent=2)+"\n").encode())
    for r in results:
        print(r["rule"],r["setting"],r["n"],r["decisions"]["5"],"gain max pp",r["largest_radius_case"]["radius_gain_pp"],"material",r["material_radius"])


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify",action="store_true")
    args=parser.parse_args()
    if args.verify:
        with tempfile.TemporaryDirectory() as d:
            build(Path(d))
            for file in OUTPUTS:
                data=(OUT/file).read_bytes()
                assert data==(Path(d)/file).read_bytes(),file
                print(hashlib.sha256(data).hexdigest(),file)
    else:
        build(OUT)
