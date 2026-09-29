# Reproducible Verification for Generalized Kerdock Codes

## Overview

This repository accompanies “The Symbol-Pair Hamming Weight Spectra of
Generalized Kerdock Codes”. It contains computational scripts and complete
machine-readable outputs for the finite cases explicitly considered in the
paper. The archived Python programs and JSON outputs are preserved unchanged.
`run_reproduction.py` provides an isolated execution interface; it does not
change the mathematical algorithms or overwrite the archived results.

## Repository structure

```text
Kerdock_Verification_GitHub/
├── README.md
├── requirements.txt
├── .gitignore
├── run_reproduction.py
├── verify.py
├── extra_checks.py
├── consistency.py
├── field_audit.py
├── review2_checks.py
├── additional_checks/
│   ├── all_primitive.py
│   └── verify_kerdock.py
└── results/
    ├── results.json
    ├── extra_results.json
    ├── consistency_results.json
    ├── all_primitive_results.json
    ├── q8_results.json
    ├── q4_m3_all_orders.json
    ├── q4_m5_rechecked.json
    ├── q8_m3_rechecked.json
    └── additional/
        └── verification-results.json
```

The runner creates `_runs/` for new computations. This release directory
contains no `_runs/`, `_audit/` or cache directories. These local-only paths
remain excluded by `.gitignore`.

## Requirements

Python 3.10 or newer and NumPy are required. The preparation checks used
Python 3.12 and NumPy 2.3.5; the tested NumPy version is pinned in
`requirements.txt`. Using the pinned dependency, Python 3.12 is recommended.
No `galois` package, API credentials or network services are used by the
computations. The remaining imports are Python standard-library or local modules.

```bash
python -m pip install -r requirements.txt
```

## Binary finite cases

For q = 2, m = 3 and m = 5, exhaustive enumeration constructs all 256 and
4,096 distinct codewords, respectively. All 32,640 and 8,386,560 unordered
codeword pairs are checked per primitive-element representative.

| m | Minimum nonzero symbol-pair Hamming weight | Minimum symbol-pair Hamming distance |
|---|---:|---:|
| 3 | 9 | 9 |
| 5 | 40 | 39 |

Thus d_SPH(K_2(4)) = 9 and d_SPH(K_2(6)) = 39. The distinction between
minimum weight and minimum distance matters because these codes are nonlinear.

Use f_3 = X^3 + X + 1 and f_5 = X^5 + X^2 + 1, with outer order
1, theta, ..., theta^(2^m-2), 0 and inner order (0,1).
`extra_checks.py` checks all Frobenius-orbit representatives (2,3 for m=3;
2,3,6,7,8,9 for m=5) and independently reconstructs the codes over Z4[X]/(f_m).
`additional_checks/all_primitive.py` separately repeats the exhaustive distance calculation for
all 6 and 30 primitive elements. Integers encode polynomial coefficients
with the constant coefficient least significant.

`results.json` also records the binary N_delta sets for m = 7,9,11, using
X^7+X+1, X^9+X^4+1 and X^11+X^2+1, respectively. These are weight and
zero-pair computations, not exhaustive pair-distance calculations for those m.

## Nonbinary finite cases

For q=4, use F4 = F2[W]/(W^2+W+1), omega = W, and
f_3 = X^3+X^2+X+omega or f_5 = X^5+X^3+X+omega.
With theta = X mod f_m, omega = theta^21 or theta^341, respectively.
The inner order is (0,1,omega,omega^2).

For q=8,m=3, use F8 = F2[W]/(W^3+W+1) and
f_3 = X^3+omega X^2+omega^5 X+omega. Here omega = theta^73 and the inner
order is (0,1,omega,omega^2,omega^3,omega^4,omega^6,omega^5).

| (q,m) | Attained weight count | Candidate weight count | Minimum nonzero weight |
|---|---:|---:|---:|
| (4,3), stated order | 35 | 40 | 204 |
| (4,5), stated order | 102 | 136 | 3264 |
| (8,3), stated order | 77 | 166 | 3640 |

For q=4,m=3, the scan considers all 36 primitive elements and both inner
orders, using omega_k = theta_k^21. Of these 72 configurations, 48 attain
35 weights and 24 attain 34; their sets differ only in weight 248.
These computations do not enumerate all nonbinary codeword pairs.

The larger weight computations use the exact normalization z=xi_0*x and
gamma=xi_1/xi_0 when xi_0 is nonzero. Every nonzero xi_0 is included through
its boundary correction; xi_0=0 is evaluated directly. This orbit reduction
still covers all q^(2m+2) parameter tuples and is not random sampling.
The programs assert the total frequency. Only the first and last inner
symbols affect these pair-weight counts: permutation blocks have one zero
and no internal zero pair. An implementation that stores only these endpoints
therefore also covers the explicitly stated inner order.

## How to reproduce the computations

Run commands from the repository root. Each command creates a new `_runs/`
directory, checks syntax and dependency imports, runs unchanged archival
algorithms and compares the generated JSON objects with archived references.
Comparisons ignore JSON formatting but include all stored values.

```bash
python run_reproduction.py smoke
python run_reproduction.py binary5
python run_reproduction.py extra
python run_reproduction.py consistency
python run_reproduction.py nonbinary
```

`smoke` checks q=2,m=3 with two implementations. `binary5` does the same for
m=5. `extra` checks binary Frobenius representatives and the independent ring
construction. `consistency` reconstructs candidate sets from archived inputs;
it does not regenerate those inputs. `nonbinary` performs the 72-order scan
and the q=4,m=5 and q=8,m=3 weight computations.

The following opt-in commands run wider or repeated exhaustive computations:

```bash
python run_reproduction.py weights
python run_reproduction.py q8
python run_reproduction.py all-primitive
python run_reproduction.py additional
```

`weights` runs all cases in `verify.py`, including binary m=11.
`q8` regenerates the original eight-ary output through `verify.py` functions.
`all-primitive` repeats binary pair enumeration for every primitive element.
`additional` runs the independent full finite-case implementation. These jobs
can require substantially more time or memory; no portable runtime is promised.

Do not invoke the archival drivers directly from the repository root:
`verify.py`, `extra_checks.py`, `consistency.py` and `additional_checks/all_primitive.py` use
relative `verification/` paths, while `review2_checks.py` and
`additional_checks/verify_kerdock.py` write beside their script files. Some drivers also start
full computations when imported. The runner preserves those conventions in
isolated copies, and checks their imports without executing the top-level loops.

## Output files

All reference JSON files below retain their original contents and names.

| File in results/ | Contents |
|---|---|
| results.json | Binary m=3,5,7,9,11 and four-ary m=3,5 weight sets, frequencies, N_delta sets and unit witnesses; binary m=3,5 distances and witnesses. Four-ary endpoint 2 and endpoint 3 are recorded separately. |
| extra_results.json | Every binary Frobenius representative, distance witnesses, minimum weights and direct ring-construction checks. |
| consistency_results.json | Candidate weight sets, unattained candidates, N_delta containment flags and theoretical Hamming-frequency totals. |
| all_primitive_results.json | Exhaustive pair-distance results and witnesses for all 6 and 30 binary primitive elements. |
| q8_results.json | Original eight-ary weight frequencies, N_delta values and unit witnesses. |
| q4_m3_all_orders.json | All 72 four-ary configurations with their attained weight and N_delta sets; this scan does not store frequency tables per configuration. |
| q4_m5_rechecked.json | Complete four-ary attained and candidate weight sets, frequencies and unattained candidates. |
| q8_m3_rechecked.json | Complete eight-ary attained and candidate weight sets, frequencies and unattained candidates. |
| additional/verification-results.json | Independent finite-case weight and N_delta computations; binary m=3,5 pair distances and word checksums. |

Fresh results and comparison reports are saved only under `_runs/`; the
reference `results/` directory is not overwritten. The preparation session
reran `smoke`, `binary5`, `extra`, `consistency` and `nonbinary` successfully.
Other archived computations were inspected and cross-compared but were not
all rerun during preparation.

## Scope and limitations

These finite computations verify only the explicitly listed finite parameter
examples under their specified field representations and coordinate orders.
They do not establish a general exact symbol-pair Hamming weight spectrum
for arbitrary q,m, or invariance of the spectrum under arbitrary coordinate
permutations. Candidate sets are supersets, not claims that every candidate
is attained. The recorded computations do not replace proofs of the general
theorems or verification of external references.
