"""Typeset the article and supplement and check their labels, bibliography and PDF links."""
from collections import Counter
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
JOBS = ("main", "low-term-supplement")
WARNINGS = re.compile(r"LaTeX Warning:|Package \S+ Warning:|(?:Over|Under)full \\[hv]box|"
                      r"Missing character:|pdfTeX warning|undefined references", re.I)


def command(args, cwd, env=None):
    run = subprocess.run(args, cwd=cwd, env=env, text=True, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, check=False)
    if run.returncode:
        raise RuntimeError(f"{' '.join(map(str, args))}\n{run.stdout}")
    return run.stdout


def build(directory):
    directory.mkdir()
    env = dict(os.environ, SOURCE_DATE_EPOCH="1795392000", FORCE_SOURCE_DATE="1",
               TZ="UTC", LC_ALL="C", TEXINPUTS=".:")
    for name in JOBS:
        # Each entry sees only its own source, the shared TeX inputs and the bibliography.
        job_dir = directory / name
        job_dir.mkdir()
        shutil.copytree(ROOT / "tex", job_dir / "tex")
        shutil.copyfile(ROOT / (name + ".tex"), job_dir / (name + ".tex"))
        command(["latexmk", "-norc", "-pdf", "-halt-on-error", "-interaction=nonstopmode",
                 "-file-line-error", name + ".tex"], job_dir, env)
        log = (job_dir / (name + ".log")).read_text(errors="replace")
        warnings = [line for line in log.splitlines() if WARNINGS.search(line)]
        if warnings:
            raise RuntimeError(f'{name} final TeX log:\n' + "\n".join(warnings))
        extracted = command(["pdftotext", name + ".pdf", "-"], job_dir)
        if "??" in extracted or "\ufffd" in extracted:
            raise RuntimeError(f"Unresolved or damaged running text in {name}")
        shutil.copyfile(job_dir / (name + ".pdf"), directory / (name + ".pdf"))
    return directory


def verify_destination(reader, key):
    dest = reader.named_destinations.get(str(key))
    if dest is None:
        raise ValueError(f"Missing named destination {key}")
    page = reader.get_destination_page_number(dest)
    if page is None or not 0 <= page < len(reader.pages):
        raise ValueError(f"Invalid destination page for {key}")


def check_labels_and_bibliography(directory, name, reader):
    job = Path(name).stem
    biblog = (directory / (job + ".blg")).read_text()
    if "Warning--" in biblog or "error message" in biblog:
        raise ValueError(f"BibTeX diagnostics in {job}\n{biblog}")
    aux = (directory / (job + ".aux")).read_text()
    cited = {key for group in re.findall(r"^\\citation\{([^}]+)\}", aux, re.M)
             for key in group.split(",")}
    entries = re.findall(r"^\\bibcite\{([^}]+)\}", aux, re.M)
    destinations = reader.named_destinations
    if not cited or cited != set(entries) or len(entries) != len(cited):
        raise ValueError(f"Mismatch between citations and bibliography entries in {job}")
    for key in cited:
        verify_destination(reader, "cite." + key)

    # Hyperref's five-field labels contain the printed page and the PDF target.
    # Other \newlabel records created by amsart hold table-of-contents widths.
    labels = [line for line in aux.splitlines()
              if line.startswith(r"\newlabel{") and not line.startswith(r"\newlabel{tocindent")]
    record = re.compile(r"\\newlabel\{([^{}]+)\}\{\{.*?\}\{(\d+)\}\{.*\}\{([^{}]+)\}\{.*\}\}")
    seen = set()
    for line in labels:
        match = record.fullmatch(line)
        if not match:
            raise ValueError(f"Unreadable cross-reference in {name}: {line}")
        key, printed, target = match.groups()
        if key in seen or target not in destinations:
            raise ValueError(f"Duplicate label or missing destination in {name}: {key}")
        seen.add(key)
        if reader.get_destination_page_number(destinations[target]) + 1 != int(printed):
            raise ValueError(f"Page of label {key} differs from its target in {name}")
    folios = [str(i) for i in range(1, len(reader.pages) + 1)]
    if reader.page_labels != folios:
        raise ValueError(f"Physical and numbered pages differ in {name}")
    for i, folio in enumerate(folios):
        target = destinations.get("page." + folio)
        if target is None or reader.get_destination_page_number(target) != i:
            raise ValueError(f"Incorrect folio destination {folio} in {name}")
    print(f"PASS {name}: {len(cited)} cited entries, {len(labels)} label pages and "
          f"{len(folios)} numbered pages agree with PDF destinations.")


def check_file_actions(reader, name):
    """Inspect reachable PDF objects, including actions outside page annotations."""
    from pypdf.generic import IndirectObject
    forbidden = {"/GoToR", "/GoToE", "/Launch", "/ImportData", "/SubmitForm", "/JavaScript"}
    pending = [reader.trailer]
    seen = set()
    while pending:
        obj = pending.pop()
        if isinstance(obj, IndirectObject):
            key = (obj.idnum, obj.generation)
            if key in seen:
                continue
            seen.add(key)
            obj = obj.get_object()
        if isinstance(obj, dict):
            if obj.get("/S") in forbidden or obj.get("/Type") == "/Filespec":
                raise ValueError(f"File or script action in {name}")
            if obj.get("/S") == "/URI" and not str(obj.get("/URI", "")).startswith("https://"):
                raise ValueError(f"Unexpected URI action in {name}")
            pending.extend(obj.values())
        elif isinstance(obj, list):
            pending.extend(obj)


def inspected_links(directory):
    """Check each PDF on its own, then check the printed references between them."""
    from pypdf import PdfReader
    palomar = json.loads((ROOT / "inputs/palomar-samples.json").read_text())["source_url"]
    for job in JOBS:
        name = job + ".pdf"
        reader = PdfReader(directory / name, strict=True)
        kinds = Counter()
        urls = set()
        check_file_actions(reader, name)
        if reader.attachments:
            raise ValueError("The computation distribution is separate from the PDFs")
        for page in reader.pages:
            for annotation in page.get("/Annots", []):
                item = annotation.get_object()
                if item.get("/Subtype") != "/Link":
                    continue
                action = item.get("/A")
                if action is None:
                    verify_destination(reader, item["/Dest"])
                    kinds["/GoTo"] += 1
                    continue
                kind = str(action["/S"])
                kinds[kind] += 1
                if kind == "/GoTo":
                    verify_destination(reader, action["/D"])
                elif kind == "/URI":
                    urls.add(str(action["/URI"]))
                    if not str(action["/URI"]).startswith("https://"):
                        raise ValueError(f"Unexpected web link in {name}")
                else:
                    raise ValueError(f"Unsupported PDF link action {kind}")
        print(f'PASS {name}: {len(reader.pages)} pages, {len(reader.named_destinations)} destinations, '
              f'{kinds["/GoTo"]} internal links and {kinds["/URI"]} HTTPS links; '
              'no file links or attachments.')
        printed = command(["pdftotext", "-raw", name, "-"], directory)
        if palomar not in urls or palomar not in "".join(printed.split()):
            raise ValueError(f"Missing printed or linked frozen Palomar source in {name}")
        lu = ("https://github.com/h-lu/classical-lienard-degree-public/tree/"
              "17a217eb2b1832f3ea75027981d93cf57b67a360/source")
        if lu not in urls:
            raise ValueError(f"Missing link to the frozen degree preprint in {name}")
        if job == "main":
            doi = "10.1007/s11425-026-2694-6"
            if "https://doi.org/" + doi not in urls or doi not in "".join(printed.split()):
                raise ValueError(f"Missing bibliographic DOI for Cen et al. in {name}")
        check_labels_and_bibliography(directory / job, name, reader)
    check_printed_locators(directory)


def check_printed_locators(directory):
    """Compare static text in one PDF with the other entry's own numbered labels after typesetting."""
    checks = {
        ("main", "low-term-supplement"): {
            "thm:lower": "Theorem {}(ii) of the electronic supplement",
            "subsec:nine-powers": "Section {} of the supplement extends its nine-term example",
            "sec:reproduce": "Section {} of the supplement",
            "subsec:full-support": "electronic supplement, Section {}.",
            "sec:trinomials": ("Section {} of the electronic supplement also proves",
                               "in Section {} of the electronic supplement. A regular value"),
        },
        ("low-term-supplement", "main"): {
            "sec:toolkit": "Article, Section {}",
            "tab:ten-comparison": "Article, Table {}",
            "subsec:other-comparisons": "Article, Section {}, stop",
            "lem:variational": "Article, Lemma {}",
            "lem:tail-inward": "Article, Lemma {}",
            "lem:hyperbolicity": "Article, Lemma {}",
            "eq:b-deriv": "Article ({})",
            "eq:b-param": "Article ({})",
            "eq:turn-derivs": "Article ({})",
            "eq:envelope-bounds": "Article ({})",
            "cor:ten-support": "Article, Corollary {}",
            "lem:compact-seed": "Article, Lemma {}",
            "lem:valley-local": "Article, Lemma {}",
            "lem:common-crest": "Article, Lemma {}",
            "prop:fold-closure": "Article, Proposition {} and Corollary",
            "cor:sign-cost": "and Corollary {}, with",
            "cor:short-lists": "and Corollary {}.",
            "prop:outer-extension": "Article, Proposition {}",
        },
    }
    for (citing, cited), locators in checks.items():
        raw = command(["pdftotext", "-raw", citing + ".pdf", "-"], directory)
        # Poppler keeps discretionary line-end hyphens, e.g. "Sec-\ntion".
        printed = " ".join(re.sub(r"(?<=\w)-\n(?=[a-z])", "", raw).split())
        aux = (directory / cited / (cited + ".aux")).read_text()
        for label, wording in locators.items():
            match = re.search(r"\\newlabel\{" + re.escape(label) + r"\}\{\{([^{}]+)\}", aux)
            wordings = (wording,) if isinstance(wording, str) else wording
            if not match or any(item.format(match[1]) not in printed for item in wordings):
                raise ValueError(f"Printed locator in {citing} for {cited} label {label}")
        total = sum(1 if isinstance(item, str) else len(item) for item in locators.values())
        print(f"PASS {citing}.pdf: {total} printed cross-document locators match {cited} numbering.")
