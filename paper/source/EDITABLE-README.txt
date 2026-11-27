Editable article, electronic supplement and computations for sparse classical Lienard systems

The standard suite contains 15 entries. Additional checks cover the fixed
degree-22 valley profile and the rational half-return data cited in the
supplementary full-support comparison. For computations alone use
source/requirements-numerical.txt; source/requirements.txt additionally
installs pypdf for the two-document typesetting checks.

The two LaTeX entries, all shared TeX files, bibliography, programs,
inputs and transcripts are in source/. The guide in source/COMPUTATION-README.txt
maps the finite calculations to the results they support. The common-crest
hypothesis and its applications are detailed in source/verification/crest-domain.md.
The primary-source locators for the return and term comparisons are in
source/verification/public-sources.md.
The printed-data check verifies the rising-profile example in Article §3.2
and all the comparison rows in their respective papers, including both
columns of lower bounds in Article Table 1. It recovers the term counts at
degrees twelve, fourteen and fifteen in Lu's construction, all three
remainders in the asymptotic term bound, and the first-order correction for
matched turns.
From the root of this extraction use the following commands.

    sha256sum -c SHA256SUMS
    python3 -m pip install -r source/requirements.txt
    scratch=$(mktemp -d)
    cd source
    python3 -B verify.py --certificates-only
    python3 -B verification/replay_certificates.py > "$scratch/all.txt"
    cmp "$scratch/all.txt" verification/certificates-full.txt
    python3 -B verification/check_stride_two.py > "$scratch/stride-two.txt"
    cmp "$scratch/stride-two.txt" verification/stride-two-fixed.txt
    python3 -B verification/check_full_support.py > "$scratch/full-support.txt"
    cmp "$scratch/full-support.txt" verification/full-support.txt
    python3 -B verification/crosscheck.py > "$scratch/ordinary.txt"
    cmp "$scratch/ordinary.txt" verification/crosscheck.txt
    python3 -B verification/check_paper_inputs.py > "$scratch/printed.txt"
    cmp "$scratch/printed.txt" verification/paper-inputs.txt
    python3 -B verification/build_papers.py --output-dir "$scratch/pdfs"

For typesetting, install pdfLaTeX, BibTeX, latexmk and Poppler. The
article uses standard amsart/AMS fonts, Latin Modern, mathtools, booktabs,
microtype, PGF/TikZ, geometry, natbib, hyperref and the abbrvnat style. Each paper
can also be built independently by running latexmk -pdf on its entry in
source/. The driver builds each entry twice in its own empty TeX tree.
References between the two documents are printed section and theorem
numbers. After typesetting, the driver checks them against each entry's
numbering and checks the PDFs for file actions. Add
--compare-dir and the directory of a reference PDF pair to verify byte
identity. The computation scripts depend only on requirements-numerical.txt.

From a clean editable extraction both ZIP files can also be reproduced with
the standard-library packer. Run in source/ with a new output directory:

    python3 -B verification/package_sources.py "$scratch/archives"

Archive membership follows the delivered source inventory and certificate
registry. The two distributions use identical numerical files; the smaller
one can be used in an environment without any TeX or PDF software.
