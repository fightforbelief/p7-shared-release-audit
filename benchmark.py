"""Exploratory single-process wall-time comparison; timing is machine dependent."""
import json
import platform
import random
from statistics import median
from time import perf_counter
from audit import Group, envelope, linear_certificate
from study import OUT, write_csv


def run():
    records = []
    for k in (8,16,32,64,128,256):
        rng = random.Random(f"P7-scaling:{k}")
        anchor = (1/k,)*k
        q = tuple(rng.uniform(.1,.9) for _ in anchor)
        groups = []
        for mass in (q,tuple(1-x for x in q)):
            rates = tuple(rng.uniform(.05,.95) for _ in anchor)
            groups.append(Group(mass,tuple(x-.03 for x in rates),tuple(x+.03 for x in rates)))
        for method in (linear_certificate,envelope):
            method(anchor,groups,.1)
            times = []
            for _ in range(5):
                start = perf_counter()
                method(anchor,groups,.1)
                times.append(perf_counter()-start)
            records.append(dict(sites=k,method=method.__name__,median_seconds=median(times),
                                repetitions=5,python=platform.python_version(),processor=platform.processor()))
    write_csv("timing.csv",records)
    print(json.dumps(records[-2:]))


if __name__=="__main__":
    run()
