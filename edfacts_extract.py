"""Extract the California public rows used in the external release audit."""
import csv
import hashlib
import io
from pathlib import Path
import sys
from zipfile import ZipFile

SOURCE_SHA256 = "b86787d850fe35a8a8200c55762ec328a7902aaa609621ac7c5b71f8e66bdf13"
EXTRACT = Path("data/edfacts_ca_math_1819.csv")
FIELDS = ("LEAID","LEANM","NCESSCH","CATEGORY","NUMVALID","PCTPROF")

def extract(path):
    if hashlib.sha256(path.read_bytes()).hexdigest()!=SOURCE_SHA256:
        raise ValueError("unexpected national archive")
    with ZipFile(path) as z, z.open("EDF_SCH_AP_MTH_1819_PUBL.csv") as raw:
        reader = csv.DictReader(io.TextIOWrapper(raw,encoding="cp1252"))
        rows = [{k:r[k] for k in FIELDS} for r in reader
                if r["FIPST"]=="06" and r["GRADE"]=="00" and r["CATEGORY"] in ("ALL","LEP")]
    if not rows:
        raise ValueError("no matching public rows")
    rows.sort(key=lambda r:(r["LEAID"],r["NCESSCH"],r["CATEGORY"]))
    with EXTRACT.open("w",encoding="utf-8",newline="") as output:
        writer = csv.DictWriter(output,fieldnames=FIELDS,lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(len(rows),hashlib.sha256(EXTRACT.read_bytes()).hexdigest())

if __name__=="__main__":
    extract(Path(sys.argv[1]))
