"""Exploratory next-year replay; all weights derived from observed applications."""
import csv
import io
import json
import zipfile
from statistics import median
from audit import envelope
from study import DATA, FILES, FIELDS, OUT, make_case, panels, write_csv

UC = ("110635","110644","110653","110662","110671","110680","110705","110714","445188")


def records(year):
    # panels() verifies the source SHA before this second read.
    with zipfile.ZipFile(DATA/f"ADM{year}.zip") as z, z.open(FILES[year][1]) as raw:
        rows = {}
        for row in csv.DictReader(io.TextIOWrapper(raw,encoding="utf-8-sig")):
            if any(row["X"+f].strip()!="R" for f in FIELDS):
                continue
            values = tuple(int(row[f]) for f in FIELDS)
            w,m,aw,am = values
            if min(w,m)>0 and 0<=aw<=w and 0<=am<=m:
                rows[row["UNITID"]] = (row["UNITID"],*values)
        return rows


def gap(weights, groups):
    rates = []
    for g in groups:
        mass = sum(b*q for b,q in zip(weights,g.share))
        rates.append(sum(b*q*r for b,q,r in zip(weights,g.share,g.low))/mass)
    return rates[0]-rates[1]


def run():
    cases,_ = panels(2022)
    panels(2023)
    old, new = records(2022), records(2023)
    cases.append(("UC",0,[old[uid] for uid in UC]))
    rows = []
    skipped = []
    for size,idx,base in cases:
        missing = [r[0] for r in base if r[0] not in new]
        if missing:
            skipped.append(dict(size=size,panel=idx,unavailable=missing))
            continue
        target = [new[r[0]] for r in base]
        b,g,_ = make_case(base)
        w,observed,_ = make_case(target)
        rho = max(0,max(1-wi/bi for bi,wi in zip(b,w)))
        actual = gap(w,observed)-gap(b,g)
        replay = gap(w,g)-gap(b,g)
        eq = envelope(b,g,rho)
        _,box,_ = make_case(base,10)
        er = envelope(b,box,rho)
        if not eq.lower.value-1e-12 <= replay <= eq.upper.value+1e-12:
            raise AssertionError("frozen replay escaped its certificate")
        rows.append(dict(size=size,panel=idx,n=len(base),rho=rho,observed=actual,
                         frozen=replay,residual=actual-replay,lo=eq.lower.value,hi=eq.upper.value,
                         round10_lo=er.lower.value,round10_hi=er.upper.value))
    write_csv("next_year.csv",rows)
    synthetic = [r for r in rows if r["size"]!="UC"]
    summary = dict(skipped=skipped,n=len(synthetic),
        outside=sum(not r["lo"]-1e-12<=r["observed"]<=r["hi"]+1e-12 for r in synthetic),
        outside_round10=sum(not r["round10_lo"]-1e-12<=r["observed"]<=r["round10_hi"]+1e-12 for r in synthetic),
        median_abs_residual_pp=median(100*abs(r["residual"]) for r in synthetic),
        median_rho=median(r["rho"] for r in synthetic),
        uc=next(r for r in rows if r["size"]=="UC"))
    (OUT/"replay.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary))


if __name__=="__main__":
    run()
