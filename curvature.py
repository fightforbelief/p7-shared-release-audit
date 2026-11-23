"""Exploratory linear-time certificate and second-order example."""
import csv
import json
from statistics import median
from audit import linear_certificate
from study import OUT, make_case, panels, write_csv


def run():
    cases = {(y,size,i):rows for y in (2022,2023) for size,i,rows in panels(y)[0]}
    with (OUT/"profiles.csv").open() as f:
        profiles = list(csv.DictReader(f))
    records = []
    for row in profiles:
        y,size,i,step,rho = int(row["year"]),row["size"],int(row["panel"]),int(row["step"]),float(row["rho"])
        b,g,_ = make_case(cases[y,size,i],step)
        lo,hi = linear_certificate(b,g,rho)
        if not lo <= float(row["lo"])+1e-12 or not hi >= float(row["hi"])-1e-12:
            raise AssertionError("curvature certificate missed")
        inner = max(0,float(row["single_hi"]))-min(0,float(row["single_lo"]))
        width = float(row["hi"])-float(row["lo"])
        records.append(dict(year=y,size=size,panel=i,step=step,rho=rho,lo=lo,hi=hi,
                     outer_excess_pp=100*(hi-lo-width),inner_shortfall_pp=100*(width-inner)))
    write_csv("curvature.csv",records)
    for rho in (.01,.1,.25,1):
        rows = [r for r in records if r["rho"]==rho and r["step"]==10]
        print(rho,median(r["outer_excess_pp"] for r in rows),max(r["inner_shortfall_pp"] for r in rows))
    return records


if __name__=="__main__":
    run()
