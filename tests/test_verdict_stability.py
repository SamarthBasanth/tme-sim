"""The panel must not tell a compound its verdict reversed when it did not.

The defect
----------
``results/verdict_inversion_census.md`` measured the shipped library on the
flagship section with the transvascular exchange term OFF vs ON, and it
measured **two different reversals** per compound:

* ``delivery_dose_verdict_inverts`` -- the panel's own headline sentence,
  delivery vs dose among the exposure-limited survivors.  **5 of 12.**
* ``dominant_limiter_inverts`` -- the four-way limiter the report prints.
  **7 of 12.**

``app/server/jobs.py`` carried ONE list: the four compounds the census's
section 4 recommends trusting.  Everything else scored was rendered as
"<compound> REVERSES its delivery/dose verdict", beside a citation of the
census.  For ``irinotecan``, ``palbociclib`` and ``ribociclib`` that is the
opposite of what was measured -- their four-way limiter flips, their
delivery/dose verdict does not -- so a reader who followed the citation found
the file contradicting the sentence that cited it.

Why a test and not a computation
--------------------------------
The lists cannot be READ from ``results/`` at runtime: ``.dockerignore``
excludes ``results/`` from the image, so a runtime read would work in the repo
and raise in the container.  They are therefore hand-carried literals -- and
this file is the only thing that stops them drifting from the run that produced
them.  Every assertion below reads
``results/verdict_inversion_census_data.json``, which is the census driver's
own output, not a second copy of the answer.

What these tests hold down
--------------------------
1. each tier is exactly what the census measured, name for name
   (:class:`TestTheTiersAreTheMeasurement`);
2. the tiers partition the scored set -- no compound in two, none in none
   (:class:`TestTheTiersPartitionTheScoredSet`);
3. the three-state verdict is reported, and the regression that motivated it
   cannot come back (:class:`TestTheVerdictIsThreeState`);
4. an unscored compound gets ``None``, never a reassuring default
   (:class:`TestAbsenceIsNotEvidence`);
5. the exposure bands quoted in the source comment are real and disjoint
   (:class:`TestTheMechanismBandsAreReal`).
"""

import json
import pathlib
import re

import pytest

from app.server.jobs import (
    VERDICT_DELIVERY_DOSE_INVERSION,
    VERDICT_LIMITER_INVERSION,
    VERDICT_LIMITER_ONLY_COMPOUNDS,
    VERDICT_REVERSING_COMPOUNDS,
    VERDICT_STABILITY_SECTION,
    VERDICT_STABLE_COMPOUNDS,
    verdict_stability,
)

CENSUS_JSON = (pathlib.Path(__file__).resolve().parents[1]
               / "results" / "verdict_inversion_census_data.json")

FINDING_JS = (pathlib.Path(__file__).resolve().parents[1]
              / "app" / "web" / "js" / "ui" / "finding.js")


@pytest.fixture(scope="module")
def prose():
    """finding.js with JS string concatenation folded away, so an assertion is
    about the sentence a reader sees and not about where the source wraps."""
    # Adjacent literals are joined whether the pair is quote-quote or
    # backtick-quote -- the tier lead is a template literal continued as a
    # plain string, and a fold that misses that seam silently stops asserting.
    return re.sub(r"[`']\s*\+\s*[`']", "", FINDING_JS.read_text())


@pytest.fixture(scope="module")
def census():
    """The census driver's own output. Absent means the assertion cannot be
    made, and a test that cannot make its assertion must not pass quietly."""
    if not CENSUS_JSON.exists():          # pragma: no cover - repo invariant
        pytest.fail(f"{CENSUS_JSON} is missing; the shipped compound tiers "
                    "cannot be checked against the run that produced them.")
    return json.loads(CENSUS_JSON.read_text())["census"]


class TestTheTiersAreTheMeasurement:
    """Each frozenset is a measured column, not a judgement call."""

    def test_reversing_is_the_delivery_dose_inverting_set(self, census):
        assert VERDICT_REVERSING_COMPOUNDS == frozenset(
            census["delivery_dose_inverting_compounds"])

    def test_limiter_only_is_limiter_inverting_minus_verdict_inverting(
            self, census):
        """The tier that did not exist before, and the reason three compounds
        were mislabelled: limiter flips, headline does not."""
        expected = (frozenset(census["dominant_limiter_inverting_compounds"])
                    - frozenset(census["delivery_dose_inverting_compounds"]))
        assert VERDICT_LIMITER_ONLY_COMPOUNDS == expected
        # Named explicitly so the regression is legible in the failure output.
        assert VERDICT_LIMITER_ONLY_COMPOUNDS == frozenset(
            {"irinotecan", "palbociclib", "ribociclib"})

    def test_stable_is_the_compounds_neither_measure_moved(self, census):
        expected = (frozenset(census["scored_compounds"])
                    - frozenset(census["dominant_limiter_inverting_compounds"])
                    - frozenset(census["delivery_dose_inverting_compounds"]))
        assert VERDICT_STABLE_COMPOUNDS == expected

    def test_the_stable_four_are_the_census_recommendation(self):
        """Section 4 names these as the only trustworthy verdicts on this
        section. The tiering must not quietly widen that list."""
        assert VERDICT_STABLE_COMPOUNDS == frozenset(
            {"vincristine", "doxorubicin", "methotrexate", "etoposide"})

    def test_both_rates_match_the_census_and_are_not_interchanged(self, census):
        """The bug in one line: these are two numbers, and the headline
        sentence goes with the smaller one."""
        assert VERDICT_DELIVERY_DOSE_INVERSION["n"] == census[
            "n_delivery_dose_inverting"]
        assert VERDICT_LIMITER_INVERSION["n"] == census[
            "n_dominant_limiter_inverting"]
        assert VERDICT_DELIVERY_DOSE_INVERSION["n"] != VERDICT_LIMITER_INVERSION["n"]
        for shipped, key in ((VERDICT_DELIVERY_DOSE_INVERSION,
                              "delivery_dose_inversion_rate"),
                             (VERDICT_LIMITER_INVERSION,
                              "dominant_limiter_inversion_rate")):
            assert shipped["of"] == census["n_scored"]
            assert shipped["percent"] == pytest.approx(
                100.0 * census[key], abs=0.05)

    def test_the_section_is_the_flagship_the_census_ran_on(self, census):
        """The census JSON records its section in prose, not as an id, so this
        pins both halves: the constant is the catalogue's flagship id, and the
        census prose describes that same section."""
        from tmesim.data.section_catalog import DEFAULT_SECTION
        assert VERDICT_STABILITY_SECTION == DEFAULT_SECTION
        config = census["config"].lower()
        for token in ("xenium", "gbm", "flagship"):
            assert token in config, (token, census["config"])


class TestTheTiersPartitionTheScoredSet:
    """A name in two tiers, or in none, is a rendering bug waiting to happen --
    'in none' renders as NOT SCORED, which reads as reassurance."""

    def test_the_tiers_are_disjoint(self):
        assert not (VERDICT_REVERSING_COMPOUNDS & VERDICT_LIMITER_ONLY_COMPOUNDS)
        assert not (VERDICT_REVERSING_COMPOUNDS & VERDICT_STABLE_COMPOUNDS)
        assert not (VERDICT_LIMITER_ONLY_COMPOUNDS & VERDICT_STABLE_COMPOUNDS)

    def test_the_union_is_exactly_what_the_census_scored(self, census):
        union = (VERDICT_REVERSING_COMPOUNDS
                 | VERDICT_LIMITER_ONLY_COMPOUNDS
                 | VERDICT_STABLE_COMPOUNDS)
        assert union == frozenset(census["scored_compounds"])
        assert len(union) == census["n_scored"] == 12

    def test_the_degenerate_compound_is_in_no_tier(self, census):
        """bevacizumab's arms are bit-identical -- its verdict is an INPUT.
        Excluded by value in the census, and it must not acquire a tier here."""
        for name in census["degenerate_zero_source"]:
            assert verdict_stability(name) is None


class TestTheVerdictIsThreeState:

    @pytest.mark.parametrize("name", sorted(VERDICT_REVERSING_COMPOUNDS))
    def test_reversing_compounds_say_reverses(self, name):
        v = verdict_stability(name)
        assert v["verdict"] == "reverses"
        assert v["stable"] is False
        assert v["delivery_dose_verdict_inverts"] is True
        assert "REVERSES" in v["note"]

    @pytest.mark.parametrize("name", sorted(VERDICT_LIMITER_ONLY_COMPOUNDS))
    def test_limiter_only_compounds_are_not_told_their_verdict_reversed(
            self, name):
        """THE REGRESSION. These three were rendered as
        '<name> REVERSES its delivery/dose verdict' beside a citation of the
        census that measured the opposite."""
        v = verdict_stability(name)
        assert v["verdict"] == "holds_with_caveat"
        assert v["delivery_dose_verdict_inverts"] is False
        assert v["dominant_limiter_inverts"] is True
        assert "REVERSES its delivery/dose verdict" not in v["note"]
        assert "KEEPS its delivery/dose verdict" in v["note"]
        # ...and it is not a clean bill either: the limiter is still named.
        assert "limiter" in v["note"]
        assert v["stable"] is False

    @pytest.mark.parametrize("name", sorted(VERDICT_STABLE_COMPOUNDS))
    def test_stable_compounds_hold_on_both_measures(self, name):
        v = verdict_stability(name)
        assert v["verdict"] == "holds"
        assert v["stable"] is True
        assert v["delivery_dose_verdict_inverts"] is False
        assert v["dominant_limiter_inverts"] is False

    def test_stable_is_true_for_the_holding_four_and_nothing_else(self):
        """``stable`` survives in stored result.json documents, so it has to
        keep meaning 'neither measure moved' rather than 'the headline
        survived'."""
        true_for = {n for n in (VERDICT_REVERSING_COMPOUNDS
                                | VERDICT_LIMITER_ONLY_COMPOUNDS
                                | VERDICT_STABLE_COMPOUNDS)
                    if verdict_stability(n)["stable"]}
        assert true_for == set(VERDICT_STABLE_COMPOUNDS)

    def test_every_scored_compound_cites_where_it_was_measured(self):
        for name in (VERDICT_REVERSING_COMPOUNDS
                     | VERDICT_LIMITER_ONLY_COMPOUNDS
                     | VERDICT_STABLE_COMPOUNDS):
            v = verdict_stability(name)
            assert v["measured_on"] == VERDICT_STABILITY_SECTION
            assert v["source"] == "results/verdict_inversion_census.md"
            assert v["compound"] == name


class TestAbsenceIsNotEvidence:

    @pytest.mark.parametrize("name", ["aspirin", "paclitaxel", "", None,
                                      "  ", "not-a-drug"])
    def test_an_unscored_compound_gets_none(self, name):
        assert verdict_stability(name) is None

    @pytest.mark.parametrize("name", ["Temozolomide", " TEMOZOLOMIDE ",
                                      "temozolomide"])
    def test_the_lookup_is_case_and_whitespace_insensitive(self, name):
        """A run whose drug is stored capitalised must not fall through to
        'not scored', which renders as reassurance."""
        assert verdict_stability(name)["verdict"] == "reverses"


class TestTheMechanismBandsAreReal:
    """The source comment claims the tiers separate on exposure retention with
    no overlap. That is the reason the rule is mechanical rather than a list of
    names, so it is measured here rather than asserted in prose."""

    @staticmethod
    def _ratios(census, names):
        return [census["compounds"][n]["mean_auc_ratio_on_over_off"]
                for n in names]

    def test_the_three_bands_do_not_overlap(self, census):
        rev = self._ratios(census, VERDICT_REVERSING_COMPOUNDS)
        mid = self._ratios(census, VERDICT_LIMITER_ONLY_COMPOUNDS)
        hold = self._ratios(census, VERDICT_STABLE_COMPOUNDS)
        assert max(rev) < min(mid), (max(rev), min(mid))
        assert max(mid) < min(hold), (max(mid), min(hold))

    def test_the_quoted_band_edges_are_the_measured_ones(self, census):
        for names, lo, hi in (
                (VERDICT_REVERSING_COMPOUNDS, 0.054, 0.140),
                (VERDICT_LIMITER_ONLY_COMPOUNDS, 0.200, 0.633),
                (VERDICT_STABLE_COMPOUNDS, 0.693, 0.927)):
            rs = self._ratios(census, names)
            assert min(rs) == pytest.approx(lo, abs=0.0005)
            assert max(rs) == pytest.approx(hi, abs=0.0005)

    def test_the_split_never_increased_exposure(self, census):
        """Check A of the census, restated where it guards the tiering: an
        inert flag would have tied the arms and made every band meaningless."""
        for name in census["scored_compounds"]:
            assert census["compounds"][name][
                "mean_auc_ratio_on_over_off"] < 1.0


class TestThePanelProseNamesTheSameTiers:
    """The rule paragraph on the panel lists the compounds by name, in three
    bands. Those names are prose in ``app/web/js/ui/finding.js`` -- there is no
    payload carrying them -- so they are the one place the tiering can drift
    back apart. Read them out of the source and compare to the frozensets,
    which :class:`TestTheTiersAreTheMeasurement` has already pinned to the
    census.
    """

    @pytest.mark.parametrize("tier,band", [
        ("VERDICT_REVERSING_COMPOUNDS", r"collapses to 0\.05 to 0\.14[^:]*:([^.]+)\."),
        ("VERDICT_LIMITER_ONLY_COMPOUNDS", r"middle band, 0\.20 to 0\.63[^:]*:([^.]+)\."),
        ("VERDICT_STABLE_COMPOUNDS", r"survives, 0\.69 to 0\.93[^:]*:([^.]+)\."),
    ])
    def test_each_band_lists_exactly_its_tier(self, prose, tier, band):
        expected = {"VERDICT_REVERSING_COMPOUNDS": VERDICT_REVERSING_COMPOUNDS,
                    "VERDICT_LIMITER_ONLY_COMPOUNDS":
                        VERDICT_LIMITER_ONLY_COMPOUNDS,
                    "VERDICT_STABLE_COMPOUNDS": VERDICT_STABLE_COMPOUNDS}[tier]
        m = re.search(band, prose)
        assert m, f"the {tier} band is no longer in the panel prose"
        named = {w.strip().strip(',').lower()
                 for w in m.group(1).replace(' and ', ', ').split(',')
                 if w.strip()}
        assert named == set(expected)

    def test_the_panel_quotes_the_delivery_dose_rate_for_its_own_sentence(
            self, prose):
        """The defect in one number. The lead sentence is about the
        delivery/dose verdict, so the rate beside it is 5 of 12 -- not the
        four-way limiter's 7 of 12, which is named separately."""
        n_dd = VERDICT_DELIVERY_DOSE_INVERSION["n"]
        n_lim = VERDICT_LIMITER_INVERSION["n"]
        of = VERDICT_DELIVERY_DOSE_INVERSION["of"]
        assert f"THIS VERDICT REVERSES FOR {n_dd} OF THE {of} SHIPPED COMPOUNDS" \
            in prose
        assert f"THIS VERDICT REVERSES FOR {n_lim} OF THE {of}" not in prose
        assert f"reverses THIS delivery/dose verdict for {n_dd} of the {of} " \
            "scored compounds" in prose
        assert f"four-way dominant limiter printed below it for {n_lim} of " \
            f"{of}" in prose

    def test_every_count_in_the_prose_is_a_tier_size(self, prose):
        """The panel spells the tier sizes out in three more places. They are
        hand-typed numbers next to hand-typed name lists, which is how the
        first one drifted."""
        of = VERDICT_DELIVERY_DOSE_INVERSION["of"]
        assert f"one of the {len(VERDICT_STABLE_COMPOUNDS)} of {of} for " \
            "which both measures hold" in prose
        assert f"for {len(VERDICT_REVERSING_COMPOUNDS)} of the {of} scored " \
            "compounds this answer is NOT DETERMINED BY THE DATA" in prose
        assert f"and for {len(VERDICT_LIMITER_ONLY_COMPOUNDS)} more the " \
            "limiter under it is not" in prose
