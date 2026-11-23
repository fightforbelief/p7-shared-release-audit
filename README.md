# Certifying changes in group rate gaps from a shared statistical release

Yunjie Guo and Hanqing Yang · Independent Researcher, Suzhou, China

[Manuscript (PDF)](paper/manuscript.pdf) · [LaTeX source](paper/manuscript.tex) · [Study design and amendments](DESIGN.md)

A disclosure-protected table can identify how much a two-group rate gap could change under different site weights. Baseline and target must use the same feasible counts. The exact algorithm returns attainable endpoints under proportional replacement of the baseline mix. The audit distinguishes a guaranteed tolerance from an attainable violation. The code includes linear-time screens, rational Bernstein bounds and an exact calculation with an attaining recipient mix.

## External validation on released school proficiency bands

The US Department of Education published all-student and English-learner mathematics proficiency bands for the 2018–19 school year. The California audit uses 3,571 schools in 417 actual districts and the native percentages, ranges and suppression markers. The estimand holds school proficiency patterns fixed. It describes possible composition effects on the learner-minus-complement gap. The information set is the ALL and LEP all-grade rows. [EDFacts source archive and documentation](https://www.eddatacenter.org/edfacts)

* At 10% replacement and a five-percentage-point tolerance, the California sharp reference certifies 192 districts and finds admissible violations in 225. The original chord screen settles 387/417, including 162 certifications. Separate-table subtraction certifies zero.
* With identical inputs and single-school witnesses, the groupwise, total-variation and reciprocal-denominator screens settle 345, 408 and 416 decisions respectively. Intersecting the three bounds settles all 417, including every case left open by the chord screen.
* The Texas protocol was fixed after the California comparison and before inspecting Texas records. Among 506 eligible districts and 2,668 schools, the chord screen settles 483, the three linear bounds settle 505, and the degree-two Bernstein bound settles all 506. The exact reference certifies 220. California and Texas share a release and school year.
* At California’s 25% replacement and 7.5-point sensitivity, Fullerton Joint Union High has a sharp lower endpoint of −6.755 points. Two recipient schools attain it. One ordinary Bernstein subdivision already certifies the threshold. The exact solver supplies the attainable endpoint and mix in a finite edge search.
* A fixed equal-school comparison gives a further limit on the role of the full count polygon. A single exploratory directional classification changes in California under nearest-percent decoding. Its margin of 0.00190 points disappears under conservative decoding.

The submission case is a short applied-statistics methods article on sharp shared-release composition bounds, attainable mixtures and the performance of standard screens. The paired results withdraw the earlier interpretation of 387/417 as an empirical advantage of the chord screen. Its second-order example characterizes the possible error from restricting recipients to a single school. Direct use of the existing total-variation inequality also achieves a second-order allowance.

[California paired comparison](results/baseline_comparison.csv) · [District results](results/baseline_districts.csv) · [Windows timings](results/baseline_timing.json) · [Texas selection and replication](results/texas_replication.json) · [Texas district results](results/texas_districts.csv) · [Fullerton subdivision and attaining mix](results/bernstein_refinement.json)

[California exact endpoints](results/edfacts_districts.csv) · [Selection and threshold comparisons](results/edfacts_summary.json) · [Original screening and sign checks](results/edfacts_screen.json) · [Original Windows verification](results/edfacts_verification.txt)

## What changes when institution totals are public?

With a 5 percentage-point tolerance, 10% replacement and rounding to tens, the sharp audit retains all 13 previously certified small panels and newly certifies 2022 small panel 5. The total is 14/20 for the exact method and 12/20 for the linear-time outer certificate. Within the original 13, the counts are 13 and 12. The other exact-only case is 2023 small panel 9. All 40 medium and large panels pass both methods. Withholding the two residual admission cells leaves these decisions unchanged. Missing group admissions are bounded by application capacity and the total.

The public-total method uses integer vertices of a count rectangle intersected with a sum band. The upper and lower endpoints are sharp under that release. Applying the older exact solver to the tightened marginal ranges gives the same threshold decisions here. The maximum endpoint difference is below 0.00064 percentage points. These panel comparisons predate the shared-table baseline and geographic validation above. The incremental decision benefit over the scan is two panels, down from four before totals were public. Rounding and panel construction remain experimental choices.

## Reproduce

The experiments ran on hanqing-win in the project-root Python 3.12.10 virtual environment. All research code uses the standard library. From a clean checkout on Windows, the commands are:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
.\.venv\Scripts\python.exe study.py
.\.venv\Scripts\python.exe replay.py
.\.venv\Scripts\python.exe curvature.py
.\.venv\Scripts\python.exe verify.py
.\.venv\Scripts\python.exe benchmark.py
.\.venv\Scripts\python.exe public_totals.py
.\.venv\Scripts\python.exe public_totals.py --verify
.\.venv\Scripts\python.exe edfacts_audit.py
.\.venv\Scripts\python.exe edfacts_audit.py --verify
.\.venv\Scripts\python.exe edfacts_screen.py
.\.venv\Scripts\python.exe edfacts_screen.py --verify
.\.venv\Scripts\python.exe edfacts_baselines.py
.\.venv\Scripts\python.exe edfacts_baselines.py --verify
.\.venv\Scripts\python.exe edfacts_baselines.py --timing
.\.venv\Scripts\python.exe texas_replication.py
.\.venv\Scripts\python.exe texas_replication.py --verify
.\.venv\Scripts\python.exe bernstein_refinement.py
.\.venv\Scripts\python.exe bernstein_refinement.py --verify
```

Each command should exit successfully before the next. On other systems the same scripts accept a Python 3.12 interpreter. The verifier checks all 1,800 attaining witnesses by direct substitution and tests interval nesting. It regenerates seven analysis files in a temporary directory, including the comparison of linear-scan and pairwise inner widths. CSV contents are checked byte-for-byte. JSON comparisons allow OS-specific line endings and environment metadata. The numerical run and independent replay used code commit `264492d`. The wall-time benchmark is machine-dependent and reported separately. To regenerate the paper figure use Python with Matplotlib and run `python paper/figures.py`. Compile the manuscript from `paper/` with `tectonic manuscript.tex`.

## Run a new audit

```python
from audit import Group, envelope, linear_certificate

anchor = (0.6, 0.4)
groups = (
    Group((0.8, 0.3), (0.60, 0.40), (0.66, 0.44)),
    Group((0.2, 0.7), (0.50, 0.32), (0.54, 0.36)),
)
exact = envelope(anchor, groups, rho=0.10)
print(exact.lower, exact.upper)
bracket = linear_certificate(anchor, groups, rho=0.10)
print(bracket.lower, bracket.inner_lower, bracket.inner_upper, bracket.upper)
```

Each group supplies its proportion at every site, then the lower and upper outcome rates. The two group proportions sum to one at each site. All endpoints use proportions; multiply by 100 for percentage points. A witness with indices i, j and fraction t replaces the designated share by a mix putting t on site i and the remainder on j. To reconstruct an attaining table choose the upper rate where the signed coefficient in Theorem 1 is positive and the lower rate elsewhere. For the lower endpoint negate the orientation. The real-arithmetic results are exact; the implementation uses double precision.

## Data and outputs

The two bundled ZIPs are the public NCES IPEDS admissions archives retrieved on 23 November 2026. The analysis uses the revised 2022 member and the sole 2023 member of the pinned ZIP. It reads the women and men application and admitted-count fields with reported-value flags. Public totals are fixed. Rounding and Laplace noise are applied by this experiment. IPEDS did not release these inputs through the simulated mechanisms. The source URL, member names and SHA-256 pins identify the supplied vintage. [NCES 2022](https://nces.ed.gov/ipeds/datacenter/data/ADM2022.zip) · [NCES 2023](https://nces.ed.gov/ipeds/datacenter/data/ADM2023.zip)

| Output | Contents |
| --- | --- |
| [Institutions](results/institutions.csv) | Selected UNITIDs, panel keys and source counts |
| [Profiles](results/profiles.csv) | All 900 cases, exact endpoints, comparison bounds and attaining pairs |
| [Private releases](results/private_releases.csv) | 600 seeded public-data simulations, true extrema, estimated bounds and coverage |
| [Primary summary](results/summary.json) | Selection flow, panel medians and simulation coverage |
| [Next-year replay](results/next_year.csv) and [summary](results/replay.json) | Actual changes, frozen-site transport, extremal ranges and excluded panels |
| [Linear scan](results/curvature.csv) | Curvature allowances expressed through outer endpoints, excess widths and rounded plug-in comparison |
| [Timing](results/timing.csv) | Five timed runs following warm-up per size and method |

The [public-total audit](public_totals.py) reads the pinned public inputs and original profile and scan files. Its [paired panel results](results/public_total_panels.csv), [summary](results/public_total_summary.json) and [Windows verification](results/public_total_verification.txt) document both releases. Output probabilities are proportions and summary widths are percentage points. The exact polygon endpoints are checked by enumerating feasible integer counts at the returned weights, and a separate three-site oracle checks the full simplex and every feasible integer table.


## EDFacts input provenance

The bundled [California public rows](data/edfacts_ca_math_1819.csv) are a column and row subset of the federal school mathematics file, snapshot 13 August 2020. The national ZIP was obtained from the Education Data Center archival mirror. To regenerate the subset, run the [extractor](edfacts_extract.py) with the downloaded ZIP path; it verifies the original SHA-256 and preserves the published labels. The extractor produced 19,360 rows. The analysis pins the resulting bytes and records every district selection loss. Passing FIPS code 48 as the second extractor argument regenerates the 15,907 [Texas public rows](data/edfacts_tx_math_1819.csv). The Fullerton case identifies the two receiving schools by their public identifiers. [National ZIP](https://storage.googleapis.com/edc-education-exports/edfacts/math-achievement-sch-sy2018-19-long.zip) · [Federal public-file documentation](https://storage.googleapis.com/edc-education-exports/edfacts/assessments-sy2018-19-public-file-documentation.docx)
