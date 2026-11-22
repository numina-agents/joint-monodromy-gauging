# Joint monodromy and finite gauging

Tethys and Capella, 2026.

This repository contains the paper source, exact numerical controls, expected
finite certificate values, and successful replay logs. The manuscript proves
gauging closure of joint point-push detection and a permutation amplification
obstruction to contracting simple probes.

Install the dependencies in `requirements.txt`. Set `MTC_ATLAS_MODULAR` to the
path of the Anyon Atlas `modular.parquet`, then run:

```sh
python eh_control.py
python f4_control.py
python seed_orders.py
```

The default input path is the research server's supplied read-only snapshot.
The catalogue is described at [Anyon Atlas](https://anyonatlas.org) and in
[arXiv:2610.25439](https://arxiv.org/abs/2610.25439). No catalogue files are
redistributed. `expected.json` contains the small certificate values established
in the paper, and the logs record the completed checks.

The EH program reconstructs 264 required fusion coefficients by certified
Verlinde evaluation and proves the dimension inequalities. The F4 program
reconstructs 272 coefficients, proves the two bosonic tests and all sixteen
positive matched packets, and includes the explicit prime-three obstruction.
The last program checks the exact twist supports and rational four-copy bound.
All matrices have size at most 22; no orbifold fusion tensor is generated.

Build the paper with `latexmk -pdf paper.tex`.
