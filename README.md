# Positive valleys and hyperbolic limit cycles in classical Liénard systems with a k-term primitive

Haibo Lu · Shanghai Institute of Technology · [ORCID](https://orcid.org/0009-0000-2717-5968)

[Main paper](paper/main.pdf) · [Electronic supplement](paper/low-term-supplement.pdf) · [Editable sources](paper/sparse-lienard-editable-source.zip) · [Numerical replay](paper/sparse-lienard-replay.zip)

## Main result

Consider the classical Liénard system

$$
\dot x=y-F(x), \qquad \dot y=-x.
$$

The paper counts the nonzero positive-integer powers in the collected primitive $F$. It does not impose a fixed degree cap. For each integer $k\ge10$, it gives a support with exactly $k$ powers and a nonempty open set of coefficients on that support such that one field has at least $4k-17$ distinct hyperbolic limit cycles. The result is a constructive lower bound.

The construction starts with a local minimum above zero in the reflected and rescaled even profile. The ordinary return and the trajectory entering from infinity sample different rising arcs. A suitable four-power block changes both return comparisons; matching the two halves at common finite launches gives sixteen additional guaranteed cycles. The paper distinguishes this term-count gain from comparisons at fixed degree, where the cost of the new exponents matters.

The electronic supplement treats small supports, including the cases $3\le k\le9$, and gives the finite construction and verification details. The article describes which conclusions follow analytically and which particular calculations are independently replayed.

## Files and verification

| File | Contents |
| --- | --- |
| [`main.pdf`](paper/main.pdf) | Main results, mechanism, proofs and comparisons |
| [`low-term-supplement.pdf`](paper/low-term-supplement.pdf) | Small-support constructions and auxiliary calculations |
| [Editable ZIP](paper/sparse-lienard-editable-source.zip) | Both LaTeX entries, shared source, bibliography, calculation inputs and instructions |
| [Replay ZIP](paper/sparse-lienard-replay.zip) | Standalone numerical and symbolic checks, with no LaTeX dependency |
| [`SHA256SUMS`](paper/SHA256SUMS) | Integrity hashes for the distributed PDFs, sources and archives |

From the repository root, verify the frozen files and run the standard suite:

~~~sh
(cd paper && sha256sum -c SHA256SUMS)
cd paper/source
python3 -B verify.py --certificates-only
~~~

Use a Python environment with the dependencies listed in [`requirements-numerical.txt`](paper/source/requirements-numerical.txt). The suite verifies 15 complete reference transcripts. For the two-document build, install the [full Python requirements](paper/source/requirements.txt), pdfLaTeX, BibTeX, latexmk and Poppler, then run the commands in the [editable build guide](paper/source/EDITABLE-README.txt). The [computation guide](paper/source/COMPUTATION-README.txt) maps checks and outputs to the results they support.

For citation or reproduction, use the frozen [`v1.0.1` snapshot](https://github.com/h-lu/sparse-lienard-limit-cycles-public/tree/v1.0.1). It contains the final manuscript and the matching supplement, source and replay materials.
