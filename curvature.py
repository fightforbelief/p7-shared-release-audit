"""Exploratory linear-time certificate and second-order example."""
import csv
from statistics import median
from audit import Group, envelope, linear_certificate
from study import OUT, make_case, panels, write_csv


def run():
    cases = {(y,size,i):rows for y in (2022,2023) for size,i,rows in panels(y)[0]}
    with (OUT/"profiles.csv").open() as f:
        profiles = list(csv.DictReader(f))
    records = []
    truth = {}
    for row in profiles:
        y,size,i,step,rho = int(row["year"]),row["size"],int(row["panel"]),int(row["step"]),float(row["rho"])
        b,g,_ = make_case(cases[y,size,i],step)
        bracket = linear_certificate(b,g,rho)
        lo,hi = bracket.lower,bracket.upper
        if step==1:
            truth[y,size,i,rho] = float(row["lo"]),float(row["hi"])
        true_lo,true_hi = truth[y,size,i,rho]
        points = []
        for group_idx,group in enumerate(g):
            rates = []
            for record in cases[y,size,i]:
                n,admitted = record[1+group_idx],record[3+group_idx]
                released = step*((2*admitted+step)//(2*step))
                rates.append(min(n,released)/n)
            points.append(Group(group.share,tuple(rates),tuple(rates)))
        plugin = envelope(b,points,rho)
        covered = plugin.lower.value<=true_lo+1e-12 and plugin.upper.value>=true_hi-1e-12
        if not lo <= float(row["lo"])+1e-12 or not hi >= float(row["hi"])-1e-12:
            raise AssertionError("curvature certificate missed")
        inner = max(0,float(row["single_hi"]))-min(0,float(row["single_lo"]))
        if abs(inner-(bracket.inner_upper-bracket.inner_lower))>1e-12:
            raise AssertionError("linear-time inner endpoints disagree")
        width = float(row["hi"])-float(row["lo"])
        records.append(dict(year=y,size=size,panel=i,step=step,rho=rho,lo=lo,hi=hi,
                     rounded_plugin_covers=int(covered),plugin_lo=plugin.lower.value,plugin_hi=plugin.upper.value,
                     outer_excess_pp=100*(hi-lo-width),inner_shortfall_pp=100*(width-inner)))
    write_csv("curvature.csv",records)
    for rho in (.01,.1,.25,1):
        rows = [r for r in records if r["rho"]==rho and r["step"]==10]
        print(rho,median(r["outer_excess_pp"] for r in rows),max(r["inner_shortfall_pp"] for r in rows))
    return records


if __name__=="__main__":
    run()
