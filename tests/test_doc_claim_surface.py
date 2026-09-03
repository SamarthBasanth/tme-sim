"""CI for the outward-facing claim surface: `README.md` and `docs/**`.

WHY THIS FILE EXISTS
--------------------
The named recurring defect in this project is CLAIM-SURFACE DRIFT: the engine
is sound and the DOCUMENTS are wrong. A claim-surface audit on 2026-08-30
enumerated 187 claim instances across the outward-facing surface and found
**148 of them in `docs/`** — a tree with, at that point, no mechanical check of
any kind. The suite was green while the surface carried 39 outright prohibited
claims, and that is not a tension: it is the finding.

Two of the three most dangerous were still live four days later:

* `docs/DEPLOYMENT.md` opened its SCOPE paragraph — the one a data-access
  committee reads to decide whether this may be deployed at all — with "Vivum
  reconstructs a tissue section and **ranks compounds against it**". That is
  the claim withdrawn on 2026-08-21, stated as what the product does.
* `docs/VERIFICATION_METHODOLOGY.md` printed two complete five-compound
  orderings side by side and cross-validated one against the other as "the
  *identical* order the simulated kill fraction produces from an independent
  computation" — a ranking presented as corroborated, on the page whose subject
  is why ranking was retired.

WHAT THIS CHECKS, AND WHAT IT DOES NOT
--------------------------------------
It checks the two patterns that are unambiguous in isolation: an explicit
ordering of two named compounds, and a sentence saying this product ranks or
compares compounds. It does NOT read intent, does not judge tone, and cannot
find a claim phrased in a way nobody has thought of. It is a floor, not the
audit — the audit found 187 things and this finds two shapes of them.

A NEGATION IS NOT AN EXEMPTION. "the ranking is withdrawn" is fine because it
contains no ordering; "the order used to be A > B" is NOT, because the reader
still reads A > B. This is the rule the project wrote for itself in
`landing-3d/lib/copy.ts`: "in any form, hedged or not, positive or negative --
including as a limit."
"""
from __future__ import annotations

import os
import re

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))

#: Every compound the engine's panels name. A name missing here is a GAP, and
#: the fix is to add it — never to loosen the patterns below.
COMPOUNDS = (
    "cisplatin", "carmustine", "lomustine", "temozolomide", "procarbazine",
    "vincristine", "vinorelbine", "doxorubicin", "methotrexate", "etoposide",
    "irinotecan", "palbociclib", "ribociclib", "bevacizumab",
)

_C = "|".join(COMPOUNDS)

#: A named compound, an ordering operator, another named compound.
#:
#: `\**` on BOTH sides and `\s` spanning newlines, both learned by this guard
#: failing its own subject test: `**cisplatin** &rarr; vincristine` slipped
#: through because the closing bold markers sit between the name and the
#: operator, and the real violation in `app/server/report.py` was
#: `<b>cisplatin &rarr;\nvincristine</b>` — wrapped across two lines, which a
#: line-at-a-time scan cannot see at all. Matches are therefore found over the
#: WHOLE FILE and mapped back to a line number afterwards.
#:
#: `->` and `&rarr;` are included because a rendered arrow is what a reader
#: sees, whatever the source spells it.
ORDERING_RE = re.compile(
    rf"\b({_C})\b\**\s*(?:>|&gt;|→|-&gt;|->|&rarr;)\s*\**\s*\b({_C})\b",
    re.IGNORECASE)

#: The product described as doing the withdrawn thing. Deliberately narrow: it
#: matches a CLAIM ABOUT THE PRODUCT, not the word "rank" in prose about why
#: ranking was retired.
PRODUCT_RANKS_RE = re.compile(
    r"\b(ranks?|ranking|orders?|compares?)\s+(the\s+)?"
    r"(compounds|candidate therapies|drugs)\b",
    re.IGNORECASE)

#: Lines that talk ABOUT the withdrawal rather than performing it. A line is
#: exempt only if it carries one of these AND no explicit ordering.
WITHDRAWAL_CONTEXT = (
    "withdraw", "withdrawn", "no longer", "used to", "retired", "prohibit",
    "is not named", "not printed", "refuses to name", "never orders",
    "does not order", "abolished", "this product will not",
    # Explicit prohibitions. `docs/architecture.md:361` reads "**do not rank
    # drugs by it**" -- an instruction NOT to do the thing, which the pattern
    # matched as an instruction TO do it. Added as whole phrases rather than
    # the bare word "not", which would exempt almost any sentence.
    "do not rank", "do not order", "do not compare", "never rank",
    "cannot rank", "must not rank",
)


#: Markdown line noise that sits between the two halves of a wrapped ordering.
#:
#: A blockquote continuation line begins `> `, and `>` is ALSO the ordering
#: operator, so a pattern permissive enough to cross it would match nonsense.
#: The text is therefore NORMALISED first — leading quote markers, list bullets
#: and table pipes stripped per line — and only then scanned.
#:
#: Found by planting a bold, blockquoted, line-wrapped ordering and watching
#: this guard stay GREEN. The unit assertion beside it used a bare newline and
#: so never exercised the case; the real violation in docs/QUICKSTART.md was
#: inside a blockquote.
_LINE_NOISE_RE = re.compile(r"^[ \t]*(?:[>|]|[-*+] |\d+\. )+[ \t]*")


def _normalise(lines):
    """Strip leading markdown structure, preserving the line count 1:1."""
    return [_LINE_NOISE_RE.sub("", ln) for ln in lines]


def _surface_files():
    out = [os.path.join(_ROOT, "README.md")]
    docs = os.path.join(_ROOT, "docs")
    for root, _dirs, files in os.walk(docs):
        for name in sorted(files):
            if name.endswith((".md", ".txt")):
                out.append(os.path.join(root, name))
    return out


def _rel(path):
    return os.path.relpath(path, _ROOT)


@pytest.fixture(scope="module")
def surface():
    docs = {}
    for path in _surface_files():
        with open(path, encoding="utf-8") as fh:
            docs[_rel(path)] = fh.read().splitlines()
    return docs


def test_the_scan_has_a_subject(surface):
    """A guard that reads nothing passes for free.

    Proves the corpus is real and that the patterns CAN fire: the compound
    vocabulary must actually occur, or every assertion below is vacuous by
    construction.
    """
    assert len(surface) >= 15, f"only {len(surface)} files scanned"
    blob = " ".join(" ".join(v) for v in surface.values()).lower()
    hits = [c for c in COMPOUNDS if c in blob]
    assert len(hits) >= 5, f"only {hits} appear anywhere — wrong corpus?"
    # And the ordering pattern must match when a real ordering is present.
    assert ORDERING_RE.search("cisplatin > doxorubicin")
    # Bold markers between the name and the operator. The first version of this
    # guard missed this and its own subject test caught it.
    assert ORDERING_RE.search("**cisplatin** &rarr; vincristine")
    assert ORDERING_RE.search("<b>cisplatin &rarr;\nvincristine</b>"), \
        "an ordering wrapped across two lines is the shape report.py had"
    # The shape that actually defeated the first version: bold, blockquoted,
    # and wrapped. Asserted through _normalise, because that is what makes it
    # visible -- testing the regex alone here would pass and prove nothing.
    quoted = ["> breaks at **cisplatin** &rarr;", "> vincristine and a spread"]
    assert ORDERING_RE.search("\n".join(_normalise(quoted))), \
        "a blockquoted wrapped ordering must be visible after normalisation"
    assert not ORDERING_RE.search(
        "\n".join(_normalise(["| dose | 5 |", "| cisplatin | 10 |"]))), \
        "a table of one compound per row is not an ordering"
    assert not ORDERING_RE.search("cisplatin was measured at 50 uM")
    assert not ORDERING_RE.search("cisplatin, doxorubicin and etoposide")


def test_no_explicit_cross_compound_ordering(surface):
    bad = []
    for rel, lines in surface.items():
        text = "\n".join(_normalise(lines))
        for m in ORDERING_RE.finditer(text):
            line_no = text.count("\n", 0, m.start()) + 1
            bad.append(f"{rel}:{line_no}: {m.group(1)} ordered against "
                       f"{m.group(2)} — {m.group(0)[:80]!r}")
    assert not bad, (
        "cross-compound ordering was withdrawn 2026-08-21 and the prohibition "
        "admits no hedge, including as a limit or as evidence:\n  "
        + "\n  ".join(bad))


def test_the_product_is_not_described_as_ranking_compounds(surface):
    bad = []
    for rel, lines in surface.items():
        for i, line in enumerate(lines, 1):
            if not PRODUCT_RANKS_RE.search(line):
                continue
            low = line.lower()
            # A line may DISCUSS the withdrawal. It may not perform it.
            if any(w in low for w in WITHDRAWAL_CONTEXT):
                continue
            bad.append(f"{rel}:{i}: {line.strip()[:110]}")
    assert not bad, (
        "the product is described as ranking or comparing compounds, with no "
        "withdrawal context on the line:\n  " + "\n  ".join(bad))


def test_a_withdrawal_sentence_is_not_treated_as_an_exemption(surface):
    """The exemption above is for CONTEXT, not for content.

    A line saying "the order used to be cisplatin > doxorubicin" contains a
    withdrawal word AND an ordering, and must still fail — the reader takes the
    ordering away regardless of the tense. Checked here rather than trusted,
    because the exemption list is the obvious way this guard could be defanged.
    """
    line = "the ranking is WITHDRAWN; it used to read cisplatin > doxorubicin"
    assert any(w in line.lower() for w in WITHDRAWAL_CONTEXT), \
        "premise: this line is exempt from the product-ranks check"
    assert ORDERING_RE.search(line), \
        "and the ordering check must still catch it — the two are independent"
