# Analysis specification, 23 November 2026

Written before inspecting IPEDS values or running the solver. This is a timestamped development specification, not a registered confirmatory study.

## Question and estimand

For an anchor application mix b, replace a share rho with any mix of the same institutions, preserving both within-site group proportions and outcome rates. Compute the attainable worst upward and downward change in the women-minus-men admission-rate gap relative to the same anchor and same latent table. Evaluate exact counts and rectangular count uncertainty. Compare against allowing anchor and target to use different tables, against one-site replacements, and against ignoring release uncertainty. No classifier training or causal inference is involved.

## Data and cases

Use NCES IPEDS fall 2022 and 2023 admissions files. Prefer revised CSV within the official ZIP. Retain reported (R) applications and admissions for both women and men, strictly positive application denominators and valid integer admission counts. Denominators refer only to these two published categories. Sort eligible units into three bins by combined application volume: [100,1000), [1000,10000), [10000,infinity). Within each bin sort by SHA-256 of the ASCII string 'P7:' followed by UNITID. Take the first 80 institutions or all if fewer; partition sequentially into panels of eight, omit incomplete panels. These synthetic portfolios test a finite-table diagnostic. Neither institutions nor panels are sampled for population inference. Repeat the same procedure in each year; cross-year overlap is expected.

## Design

Use observed application shares for the anchor. Primary profile rho = 0.10, with 0, 0.01, 0.25 and 1 as boundary/sensitivity checks. Use an exact-count control and a hypothetical release rounding admitted counts to the nearest multiple of 10 or 100 (half ties upward, intersect with feasible integers). Applications stay public. Interval endpoints quantify ambiguity conditional on the release. Rounding has no differential privacy guarantee.

Also emulate a public-denominator outcome-private Laplace release at epsilon 1 on the first panel in each size bin and year. Use 100 independent seeded public-data simulations per panel at rho=0.10. Compare an honest simultaneous 95% noise-inversion box with clipped plug-in counts and independently ranged anchor/target. Report frequencies against exact-table extrema, uncertainty width and Monte Carlo precision. The seeded simulation is not a production privacy mechanism.

Report panel medians, ranges or frequencies with explicit denominators. Include an independent exhaustive rate-corner and grid/analytic mixture oracle on small fixtures, interior and degenerate examples and attainability checks. Disclose failed comparisons, differences from this specification and selection limits.

## Exploratory additions

The centered total-variation upper bound was added before running the study as a comparator that accounts for cancellation. After the primary run, we added a next-year diagnostic. Fix the 2022 panels, retain complete panels with valid 2023 records for every member, and add the nine University of California undergraduate campuses as a named system. Compare the observed 2023 change with the 2022 frozen-site projection at the observed 2023 weights and with its envelope at the smallest admissible replacement fraction. Neither the temporal replay nor the named case is confirmatory. Cross-year results may reflect changes in group shares, outcomes or reporting. The UC case uses application counts per campus.

The small primary-radius error of single-site replacement suggested a second exploratory analysis. We derived a uniform curvature allowance for the single-site audit and evaluated its enclosure and excess width on all existing profiles. This analysis used the existing data without redefining the primary comparisons.

## Public institution totals (follow-up specified before computing endpoints)

Reuse all 60 original panels, their baseline weights, the two-group estimand, 10% replacement, half-up rounding to tens and the inclusive ±5 percentage-point decision. Primary comparison: publish exact institution-wide application and admission totals, public application denominators, and rounded admissions for every numerically available group (women, men, another gender and gender unknown). Preserve missing cells. A missing admissions value carries only the feasibility restrictions from applications and the total. The combined other-group application capacity limits their combined admissions. Do not reconstruct a two-group total by subtracting the original unrounded residual admission counts. A second release withholds both residual admission cells entirely while preserving the total and the combined residual capacity.

Pair each decision with its original year/size/panel key. Report retention among the original 13 certified small panels, newly certified panels, and exact versus coupled linear-time certificates both within those 13 and in all 20 small panels. Nesting predicts that all 13 must survive; the decision increments are the empirical question. Evaluate medium and large panels and report interval widths as secondary outcomes. Use identical site constraints for both algorithms. Validate projected integer feasibility, compare small cases to exhaustive table enumeration, check the sharp endpoints against witnesses and the unconstrained baseline, and reproduce the new outputs. This is a post-primary paired analysis.
