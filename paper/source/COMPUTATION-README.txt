Finite computations for Positive valleys and hyperbolic limit cycles in classical Liénard systems with a k-term primitive

Use this directory by itself, the root of the standalone computation ZIP,
or source/ in the editable distribution. Python, mpmath, python-flint and
SymPy suffice. The standalone pins are in requirements.txt; their editable
counterpart is requirements-numerical.txt. From a standalone extraction:

    sha256sum -c SHA256SUMS
    python3 -m pip install -r requirements.txt
    python3 -B verify.py --certificates-only
    python3 -B code/positive_valleys.py

The additional fixed degree-22 profile has its own interval check:

    scratch=$(mktemp -d)
    python3 -B verification/check_stride_two.py > "$scratch/stride-two.txt"
    cmp "$scratch/stride-two.txt" verification/stride-two-fixed.txt

The direct positive-valleys run has a standalone reference in verification/positive-valleys-direct.txt.

The standard suite runs all 15 entries in inputs/certificates.json, including
the finite profile, refined quadrature and printed seed checks. Each complete
transcript (command, stdout, exit status) must match checks/NAME.txt.
Keep assertions enabled. The direct positive-valleys command repeats entry 10.
For the additional floating-point comparison run:

    scratch=$(mktemp -d)
    python3 -B verification/crosscheck.py > "$scratch/ordinary.txt"
    cmp "$scratch/ordinary.txt" verification/crosscheck.txt

Map of statements to executable inputs and reference outputs
-----------------------------------------------------------
In the list below, prefix each numbered command with python3 -B.
"Supplement" means the electronic supplement to this article. The named
output in the right column is checks/NAME.txt. The auxiliary ledger
inputs/supports.json supplies the exponents, comparison and target data.

Article or supplement location       Command / reference name

 1  Supplement Lemma 3.1, Table 2    code/trinomial.py
    Exact rational gap enclosures    trinomial

 2  Supplement §4.1                  code/four.py
    Three seed and three outer       four-term
    action intervals

 3  Article Lemma 3.1(iii),          code/balanced_seed.py --degree 10
    Corollary 5.2; supplement §4     degree-ten-seed
    Seed actions at 20,100,200. For the peak H_b<14, the isolated radius
    and the curvatures used in the two 16-cycle examples, also run entry 9.

 4  Supplement §6.2                  code/crest.py --new-odd 1
    Five-term walls of degrees       five-crest
    16384 and 20000000; dyadic bisection and 95-decimal-digit interval arithmetic

 5  Supplement §6.2                  code/crest.py --new-odd 3
    Exact rational three-root        seven-crest
    block before the eventual seven-term crest

 6  Article Lemma 3.1(iii),          code/distant_fold.py
    supplement §§5.1–5.2             six-fold
    Rational bounds for seed curvatures and the degree-eighteen cancellation bracket

 7  Supplement §6.3                  code/two_profiles.py --terms 8
    Seven rising-node bounds         eight-profile
    and wall estimates at degree 2^38

 8  Article §§3.2, 4.2,             code/bounded_outer.py
    supplement §6.4                 nine-bounded
    Exact cubic margins and the seven normalized floors in the sharpness
    witness, plus the finite degree-eighteen fold comparisons

 9  Article Lemma 3.1(iii);          code/baseline_checks.py
    supplement §§3–5                 baseline-rational
    Rational gap bounds, seed radius, peak and curvatures, three seed and
    three outer actions, terminal cancellation. This uses a separate
    interval implementation for the overlapping checks in entries 1–3,6.

10  Article Lemmas 3.1,3.4;          code/positive_valleys.py
    §4.3, Tables 2–3                 positive-valleys
    Rational prefix tests and support accounting, 768-bit Arb enclosures
    of the high-degree seed, critical points (including 7/8 < z_p), two
    seven-sign positive lists, and five depth tests. The 106498-cell
    improper-integral bounds include both analytic tails. All-prefix
    signs in the general theorem also have a paper proof.

11  Article Theorem 1.1, Table 1;    code/check_ledger.py
    supplement Theorem 1.1, §6.4    support-ledger
    Checks the seven low-term support sizes, nonzero coefficients of their
    radial averages, seven odd extensions and seven rational cubic margins.
    For degree 67088 it checks the exponents and the target signs in
    inputs/supports.json. It checks all four residue counts for 1 to 32 blocks.

12  Article Proposition 3.2,         code/return_identities.py
    Lemma 4.2, Appendix A;           return-identities
    supplement §2.1
    Symbolic elimination identity, moving-endpoint derivatives on a
    polynomial example, averaged binomial drift and exterior-oval derivative.

13  Article §4.2 and Corollary 5.2   verification/check_mechanism.py
    mechanism; stdout also in verification/mechanism-replay.txt
  512-bit Arb. Article §4.2 and Corollary 5.2: finite matrix and affine
  descent tests; two positive lists and three improper integrals on the
  degree-26 profile (65536 cells).

14  Article §4.3 and Corollary 5.2   verification/check_quadrature.py
    refined-quadrature; stdout also in verification/refined-quadrature.txt
  Independent whole-cell integration and tail estimates for all eight
  depth intervals. 212996 cells at 768 bits and 131072 cells at 512 bits.

15  Article Lemma 3.1(iii)           code/seed_table.py
    printed-seed. Exact factorization and antiderivative identities, endpoint
    inequalities and interval evaluation from inputs/degree-ten-seed.json.

Additional interval check outside the 15-entry suite
----------------------------------------------------
verification/check_stride_two.py -> verification/stride-two-fixed.txt
  At d=10, p=2, eta=2^-20 and q=2^-10, encloses both positive-height
  lists and J at three depths. Reciprocal substitution maps the whole
  improper integral to [0,1]. Each enclosure uses 131072 whole cells
  and 512-bit Arb. Article Corollary 5.2 gives the low-degree cycle count.

Full-support literature comparison in Supplement Section 1.1
------------------------------------------------------------
verification/check_full_support.py -> verification/full-support.txt
  The input inputs/palomar-samples.json gives the pinned unsigned preprint,
  SHA-256 of its TeX source and its outward-rounded Appendix A data. This
  check uses exact fractions to compare the fourteen printed intervals and
  check coefficient nonvanishing and the order of the seven launch heights.
  The full-support counts use the preprint's persistence and outer-barrier
  results. Hyperbolicity at nearby coefficients uses Article Lemma 2.3.
  The degree-13, ten-term row uses the four odd extensions of Article
  Proposition A.1. With --source-tex PATH the program also checks the
  source digest, all six coefficients, seven heights and 28 interval
  endpoints against Appendix A of the frozen TeX file.
  The smallest separation of the printed intervals is 2300260502/10^42.
  The original validated ODE integrations belong to the cited preprint.
  The public TeX is pinned to commit
  7da225327c6cd62a3dc63926525377764ea8f454, file paper/main.tex.

    python3 -B verification/check_full_support.py

Diagnostic outside the certificate suite
----------------------------------------
verification/crosscheck.py -> verification/crosscheck.txt
  Exact support arithmetic and ordinary 65-decimal-digit mpmath evaluation
  of four determinants, the threshold constant and three integrals.
  These floating-point values provide a numerical comparison.

The analytic estimates and matching prove Theorem 1.1, including
at least 26 cycles on ten terms for sufficiently large odd n0. Entry 10
checks the fixed seed and valley inequalities for n0=257 in Corollary 4.6,
whose prescribed support has degree 67088. For the two sixteen-cycle
examples, the valley signs follow by asymptotics in both degrees 22 and 26.
They use the same rational degree-ten seed inputs. At the same dyadic parameters, separate Arb runs enclose both profiles,
entry 13 for degree 26 and check_stride_two.py for degree 22.
The proofs choose the dilations, the common matching coefficients and the
final regular perturbation. The program outputs enclose the seed actions
and local profile integrals. Section 4.2 proves the passage to actual
finite orbits and Section 4.3 counts their disjoint brackets and marks.
The editable distribution's verification/check_paper_inputs.py also checks
the article's ten-term table and the supplement's seven small-support rows.
It independently integrates the rising-profile variation in Article §3.2
and compares the seven exact integer parts with the printed table. It evaluates the exponential
bases at the three subunit samples of Lemma 3.4 and at w=7/10.
At that additional depth the squared base is 15/14 while H/m grows as
(21/20)^n. The three original samples have squared bases 81/128, 9/8,
and 6/7, respectively.
Published full-support bounds (LMP, LL, DMH) and the preprint's exact-degree
Theorem 1.1 and Supplement Proposition S1.2 (for k >= 6) supply the
full-support comparisons. The separate degree-12 column follows from the
preclosure ordinary tests in §§5.2–5.3 of the preprint. The printed-data
check computes the one-radius moment, the two tail signs, and both
preclosure bracket counts. For the degree-14 and degree-15 rows of
Proposition 5.1 it counts 12 and 13 allowed terms and 23 and 24 cycles,
respectively. It also checks the finite floor bound obtained by truncation,
the printed critical-point and positive
enclosures in Article §4.3, and checks the turning-height correction
against (S5.12)–(S5.13) by eliminating the first-order shift.
