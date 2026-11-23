"""Frozen Texas replication of the primary EDFacts screen comparison."""
import argparse
import csv
import hashlib
import io
import json
from collections import Counter

from audit import linear_certificate
from edfacts_audit import audit, load
from edfacts_baselines import EPS, NAMES, OUT, decision, inputs, radius
from edfacts_extract import TEXAS
from shared_baselines import bounds

PIN = "04d1c82d85bfbfcf063bd5e4ec10dc2b13873dc977b817df34eaa99328a3dc90"


def build():
    selected,flow = load(False,TEXAS,PIN)
    buf = io.StringIO(newline="")
    fields = ("district","schools","method","inner_lo","inner_hi","outer_lo","outer_hi",
              "joint_lo","joint_hi","proj_lo","proj_hi","decision")
    writer = csv.DictWriter(buf,fieldnames=fields,lineterminator="\n")
    writer.writeheader()
    counts = {m:Counter() for m in NAMES}
    exact_pass=0
    projection_flips=0
    for district,(_,sites) in sorted(selected.items()):
        anchor,groups=inputs(sites)
        jl,ju,pl,pu=audit(sites,"0.10")
        passed=radius((jl,ju))<=.05+EPS
        exact_pass+=passed
        projection_flips+=passed!=(radius((pl,pu))<=.05+EPS)
        vals=bounds(anchor,groups,.1)
        chord=linear_certificate(anchor,groups,.1)
        vals["chord"]=(chord.lower,chord.upper)
        il,iu=vals["inner"]
        assert jl<=il+EPS and iu<=ju+EPS
        assert abs(il-chord.inner_lower)<EPS and abs(iu-chord.inner_upper)<EPS
        old=decision((il,iu),vals["chord"],.05)
        for m in NAMES:
            lo,hi=vals[m]
            d=decision((il,iu),(lo,hi),.05)
            assert lo<=pl+EPS and pu<=hi+EPS and pl<=jl+EPS and ju<=pu+EPS
            assert d=="open" or (d=="pass")==passed
            counts[m][d]+=1
            counts[m]["chord_only"]+=old!="open" and d=="open"
            counts[m]["comparator_only"]+=d!="open" and old=="open"
            writer.writerow(dict(zip(fields,(district,len(sites),m,il,iu,lo,hi,jl,ju,pl,pu,d))))
    summary=dict(flow=flow,exact_pass=exact_pass,projection_flips=projection_flips,
                 methods={m:{k:c[k] for k in ("pass","violation","open","chord_only","comparator_only")}
                          for m,c in counts.items()})
    return {"texas_districts.csv":buf.getvalue(),"texas_replication.json":json.dumps(summary,indent=2)+"\n"}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify",action="store_true")
    args=parser.parse_args()
    for name,s in build().items():
        data=s.encode()
        if args.verify:
            assert (OUT/name).read_bytes()==data,name
        else:
            (OUT/name).write_bytes(data)
        print(hashlib.sha256(data).hexdigest(),name)
        if name.endswith(".json"):
            print(s)
