"""
Black-box tests for parse_solvent.py, based solely on observing real data
(the `solvent` / `solvent_clean` columns of csd_all.csv).

These tests encode the DESIRED behavior, not necessarily the script's
current behavior. Some are expected to fail (xfail) until the identified
bugs are fixed.
"""
import pytest

from parse_solvent import parse_value, _UNKNOWN


# ---------------------------------------------------------------------------
# 1. Simple solvents: full name, abbreviation, varied casing
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("raw,expected", [
    ("ethanol", "ethanol"),
    ("methanol", "methanol"),
    ("Ethanol", "ethanol"),
    ("ETHANOL", "ethanol"),
    ("water", "water"),
    ("hexane", "hexane"),
    ("toluene", "toluene"),
    ("acetone", "acetone"),
    ("chloroform", "chloroform"),
    ("dichloromethane", "dichloromethane"),
    ("tetrahydrofuran", "tetrahydrofuran"),
    ("benzene", "benzene"),
    ("pentane", "pentane"),
    ("diethyl ether", "diethyl ether"),
])
def test_simple_solvent_name(raw, expected):
    assert parse_value(raw) == expected


@pytest.mark.parametrize("raw,expected", [
    ("EtOH", "ethanol"),
    ("MeOH", "methanol"),
    ("THF", "tetrahydrofuran"),
    ("DCM", "dichloromethane"),
    ("Et2O", "diethyl ether"),
    ("MeCN", "acetonitrile"),
    ("DMSO", "dimethyl sulfoxide"),
    ("DMF", "dimethylformamide"),
])
def test_abbreviations_expanded(raw, expected):
    assert parse_value(raw) == expected


@pytest.mark.parametrize("raw,expected", [
    ("diethylether", "diethyl ether"),
    ("n-hexane", "hexane"),
    ("hexanes", "hexane"),
    ("di-n-butyl ether", "dibutyl ether"),
    ("light petroleum", "petroleum ether"),
])
def test_synonyms_normalized(raw, expected):
    assert parse_value(raw) == expected


# ---------------------------------------------------------------------------
# 2. Mixtures: varied separators, alphabetical sort, deduplication
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("raw,expected", [
    ("dichloromethane/hexane", "dichloromethane/hexane"),
    ("hexane/dichloromethane", "dichloromethane/hexane"),
    ("hexane-toluene", "hexane/toluene"),
    ("water / dichloromethane", "dichloromethane/water"),
    ("hexane/methanol/dichloromethane", "dichloromethane/hexane/methanol"),
    ("methanol/water/diethyl ether/tetrahydrofuran",
     "diethyl ether/methanol/tetrahydrofuran/water"),
])
def test_mixture_sorted_alphabetically(raw, expected):
    assert parse_value(raw) == expected


def test_mixture_deduplicates_identical_solvents():
    assert parse_value("ethanol/ethanol") == "ethanol"
    assert parse_value("ethanol/Ethanol") == "ethanol"


@pytest.mark.parametrize("raw,expected", [
    ("acetone/heptane(1:1)", "acetone/heptane"),
    ("acetone/water 1:1", "acetone/water"),
    ("chloroform/benzene (4:1)/diethyl ether",
     "benzene/chloroform/diethyl ether"),
])
def test_mixture_ratio_stripped(raw, expected):
    assert parse_value(raw) == expected


# ---------------------------------------------------------------------------
# 3. Numeric/chemical noise: percentages, molarity, temperature
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("raw,expected", [
    ("95% ethanol", "ethanol"),
    ("96% ethanol", "ethanol"),
    ("80% ethanol/water", "ethanol/water"),
    ("95% ethanol/ethyl acetate", "ethanol/ethyl acetate"),
])
def test_percentage_stripped(raw, expected):
    assert parse_value(raw) == expected


@pytest.mark.parametrize("raw,expected", [
    ("pentane/toluene at -35 deg.C", "pentane/toluene"),
    ("Crystals grown from Et2O-hexane at 293 K", "diethyl ether/hexane"),
    ("Crystals grown from THF-hexane at 293 K", "hexane/tetrahydrofuran"),
])
def test_temperature_mentions_stripped(raw, expected):
    assert parse_value(raw) == expected


# bug found at runtime: the "cooling X to -Y deg.C" construction breaks
# extraction (-> None) whereas "X at -Y deg.C" works on similar entries
# just above.
@pytest.mark.xfail(reason="known bug: 'cooling X to -Y deg.C' construction not parsed")
def test_cooling_to_temperature_construction_stripped():
    assert parse_value("cooling dichloromethane to -40 deg.C") == "dichloromethane"


# known bug: a valid solvent mixed into a "X in Y" description isn't
# recognized. This encodes the desired behavior (both components
# extracted); the test should fail until the script is fixed.
@pytest.mark.xfail(reason="known bug: '<solvent> in <solvent>' construction not parsed")
@pytest.mark.parametrize("raw,expected", [
    ("40% ethyl acetate in hexanes", "ethyl acetate/hexane"),
    ("5% ethyl acetate in hexanes ", "ethyl acetate/hexane"),
])
def test_percentage_in_solvent_construction(raw, expected):
    assert parse_value(raw) == expected


# ---------------------------------------------------------------------------
# 4. Acids: normalized to full name, molarity stripped
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("raw,expected", [
    ("0.1M hydrochloric acid", "hydrochloric acid"),
    ("6 M hydrochloric acid", "hydrochloric acid"),
    ("0.1M HCl", "hydrochloric acid"),
    (">12 M hydrochloric acid", "hydrochloric acid"),
    ("6 M HBr", "hydrobromic acid"),
    ("2M sulfuric acid", "sulfuric acid"),
    ("0.2M nitric acid", "nitric acid"),
    ("7% acetic acid", "acetic acid"),
])
def test_acid_normalized(raw, expected):
    assert parse_value(raw) == expected


@pytest.mark.xfail(reason="known bug: 'aqueous X' drops the water instead of keeping it as a mixture component")
def test_aqueous_acid_keeps_water_component():
    assert parse_value("aqueous HCl") == "hydrochloric acid/water"


@pytest.mark.parametrize("raw,expected", [
    ("1M hydrochloric acid/methanol", "hydrochloric acid/methanol"),
    ("6M hydrochloric acid/N,N-dimethylformamide ",
     "dimethylformamide/hydrochloric acid"),
    ("hydrobromic acid/ethanol", "ethanol/hydrobromic acid"),
])
def test_acid_in_mixture(raw, expected):
    assert parse_value(raw) == expected


# ---------------------------------------------------------------------------
# 5. Aqueous markers: "aqueous" alone -> water (already correct), but
#    "aqueous X" should keep the water as a mixture component (known bug:
#    the script drops the water instead of keeping it).
# ---------------------------------------------------------------------------

def test_aqueous_alone_means_water():
    assert parse_value("aqueous") == "water"


@pytest.mark.xfail(reason="known bug: 'aqueous X' drops the water instead of keeping it as a mixture component")
@pytest.mark.parametrize("raw,expected", [
    ("aqueous ethanol", "ethanol/water"),
    ("aqueous acetone", "acetone/water"),
    ("aqueous/chloroform/methanol", "chloroform/methanol/water"),
])
def test_aqueous_marker_keeps_water_component(raw, expected):
    assert parse_value(raw) == expected


# ---------------------------------------------------------------------------
# 6. Deuterated solvents: isotopic information must be kept, with a
#    standard and consistent NMR notation.
# ---------------------------------------------------------------------------

@pytest.mark.xfail(reason="known bug: deuteration is dropped instead of being normalized")
@pytest.mark.parametrize("raw,expected", [
    ("deuterated benzene", "benzene-d6"),
    ("deuterated chloroform", "chloroform-d"),
    ("deuterated methanol", "methanol-d4"),
    ("deuterated dichloromethane", "dichloromethane-d2"),
])
def test_deuterated_solvent_keeps_isotope_info(raw, expected):
    assert parse_value(raw) == expected


def test_deuterated_water_already_correct():
    assert parse_value("dideuterated water") == "water-d2"


def test_deuterated_explicit_suffix_preserved():
    assert parse_value("benzene-d6") == "benzene-d6"


@pytest.mark.xfail(reason="known bug: the abbreviated '-d6' form is not normalized consistently everywhere")
def test_deuterated_mixture_consistent_suffix():
    # benzene-d6 and deuterated toluene in the same mixture: the isotopic
    # notation must be preserved for every affected component.
    assert parse_value("toluene/ deuterated benzene") == "benzene-d6/toluene"


@pytest.mark.xfail(reason="known bug: 'ammonia solution or deuterated water' loses the deuteration info")
def test_deuterated_water_in_alternative_phrasing():
    assert parse_value("ammonia solution or deuterated water") == "ammonia/water-d2"


# ---------------------------------------------------------------------------
# 7. Descriptive sentences / non-solvent text -> None
#    (only when NO solvent name appears in the text)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("raw", [
    "solution layering",
    "slow evaporation from filtrate",
])
def test_pure_descriptive_text_returns_none(raw):
    assert parse_value(raw) is None


# known bug: these sentences nonetheless contain explicit solvent names
# (pentane, dichloromethane, methanol, water, N-methyl-2-pyrrolidone...)
# buried in descriptive text. The script returns None instead of extracting
# the mentioned solvent(s), unlike other textual-noise cases (temperature,
# percentage) that are correctly ignored elsewhere.
@pytest.mark.xfail(reason="known bug: descriptive text around a valid solvent name blocks extraction")
@pytest.mark.parametrize("raw,expected", [
    ("layering of pentanes on a concentrated dichloromethane solution of the complex",
     "dichloromethane/pentane"),
    ("Layered in N-methyl-2-pyrrolidone", "n-methyl-2-pyrrolidone"),
    ("methanol, then crystals immersed in water", "methanol/water"),
    ("Crystals were obtained at room temperature by slow diffusion of pentane "
     "into a concentrated solution of the compound dissolved in dichloromethane.",
     "dichloromethane/pentane"),
    ("hexane/sublimation", "hexane"),
])
def test_descriptive_text_with_named_solvent_still_extracts_it(raw, expected):
    assert parse_value(raw) == expected


# same bug, variant with the "aqueous" marker (see section 5): the implicit
# "water" solvent should be extracted instead of returning None.
@pytest.mark.xfail(reason="known bug: 'aqueous' in descriptive text doesn't produce 'water'")
@pytest.mark.parametrize("raw,expected", [
    ("evaporation from aqueous buffer solution", "water"),
    ("evaporation from aqueous solution", "water"),
    ("High-pressure crystallisation from aqueous solution", "water"),
])
def test_descriptive_text_with_aqueous_marker_extracts_water(raw, expected):
    assert parse_value(raw) == expected


@pytest.mark.parametrize("raw,expected", [
    ("Re-crystallisation from solvent: dimethylsulfoxide", "dimethyl sulfoxide"),
    ("Re-crystallisation from solvent:  acetonitirle", "acetonitrile"),
    ("Re-crystallisation from solvent: THF", "tetrahydrofuran"),
    ("Re-crystallisation from solvent: MeOH, DCM", "dichloromethane/methanol"),
])
def test_recrystallisation_prefix_does_not_block_extraction(raw, expected):
    assert parse_value(raw) == expected


# known bug found while reviewing _resolve(): the "full name (own abbreviation)"
# unwrap (e.g. "n,n-diethylformamide (def)" -> "diethylformamide") is applied
# to the component BEFORE the leading from/of/containing prefix strip, instead
# of after. When the prefix is stripped at the whole-string level (e.g. by
# "Re-crystallisation from solvent:") this doesn't matter, but when a "from "
# ends up attached to a component after a mixture split, the abbreviation
# never gets unwrapped and the component stays unresolved.
@pytest.mark.xfail(reason="known bug: '(abbreviation)' unwrap doesn't apply after a component-level prefix strip")
def test_abbreviation_unwrap_applies_after_component_prefix_strip():
    assert parse_value("hexane/from n,n-diethylformamide (def)") == "diethylformamide/hexane"


# ---------------------------------------------------------------------------
# 8. Missing / undetermined values
#    - blank / whitespace / None / "?" -> no info -> None
#    - explicit "undetermined" markers -> _UNKNOWN sentinel
#      (triggers the script's fail-fast mode, by design)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("raw", ["", " ", None, "?"])
def test_blank_value_returns_none(raw):
    assert parse_value(raw) is None


@pytest.mark.parametrize("raw", ["n/a", "N/A", "not determined"])
def test_explicit_not_determined_returns_unknown_sentinel(raw):
    assert parse_value(raw) is _UNKNOWN


# ---------------------------------------------------------------------------
# 9. Robustness: spaces, tabs, stray punctuation
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("raw,expected", [
    ("dichloromethane/n-hexane solvate\t", "dichloromethane/hexane"),
    ("  ethanol  ", "ethanol"),
    ("ethanol\n", "ethanol"),
    ("hexane/ benzene/ dichloromethane", "benzene/dichloromethane/hexane"),
])
def test_whitespace_and_stray_characters_ignored(raw, expected):
    assert parse_value(raw) == expected


# ---------------------------------------------------------------------------
# 10. Regressions fixed on 2026-06-16, following the resolution of the
#     "Re-crystallisation from solvent:" bug (see project_parse_solvent
#     memory): this uncovered real cases previously blocked upstream by a
#     substring false positive on NON_SOLVENT_KW/DESCRIPTION_KW.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("raw", [
    "Slow Evaporation Solution Growth (SEST)",
    "Solvent exchange",
    "in situ on the diffractometer at a temperature of 216 K",
    "slow evaporation solution technique",
    "slow evaporation within MOF",
    "obtained from an aqueous solution containing L-asparagine and thioacetamide",
    "obtained from an aqueous solution containing 1.5 mol% L-proline",
    "obtained from an aqueous solution containing ammonium bromide",
    "Slow Evaporation/Photolytic",
    "Slow evaporation of the liquor mother.",
])
def test_pure_method_description_returns_none_regression(raw):
    assert parse_value(raw) is None


@pytest.mark.parametrize("raw,expected", [
    ("grown from slow evaporation of a hexane/ethyl actetate solution", "ethyl acetate/hexane"),
    ("evaporation of acetonitrile (RT)", "acetonitrile"),
    ("slow evaporation from dicloromethane", "dichloromethane"),
    ("Slow evaporation from a concentration tetrahydrofuran solution", "tetrahydrofuran"),
    ("Slow evaporation of pentante", "pentane"),
    ("slow evaporation from nitromethane/methanol (v:v 4:1)", "methanol/nitromethane"),
    ("slow evaporation of solvent (acetonitrile)", "acetonitrile"),
    ("Crystallized from slow evaporation of soluiton of dichloromethane", "dichloromethane"),
    ("crystallised by slow cooling from pentanol solution", "1-pentanol"),
    ("Slow evaporation from methanol (to dryness).", "methanol"),
    ("Slow evaporation of mix solvent of DCM and EtOH containg nanocluster 3", "dichloromethane/ethanol"),
    ("Slow evaporation of mds088 from a mixture of methanol and benzonitrile.", "benzonitrile/methanol"),
    ("Grown in the excess p-xylene", "p-xylene"),
    ("slow evaporation out of ethyl acetate/ hexane", "ethyl acetate/hexane"),
])
def test_method_description_noise_stripped_around_valid_solvent(raw, expected):
    assert parse_value(raw) == expected


@pytest.mark.parametrize("raw,expected", [
    ("Re-crystallisation from solvent: 1,2-dibromoethane", "1,2-dibromoethane"),
    ("Re-crystallisation from solvent: N,N-diethylformamide (DEF)", "diethylformamide"),
    ("Re-crystallisation from solvent: N-formylpiperidine (NFP)", "n-formylpiperidine"),
    ("Re-crystallisation from solvent: ethanole", "ethanol"),
    ("Re-crystallisation from solvent: Tol", "toluene"),
    ("Re-crystallisation from solvent: CHCl3_and_Hx", "chloroform/hexane"),
    ("Re-crystallisation from solvent: dichloromethane aand acetonitrile mixture", "acetonitrile/dichloromethane"),
    ("Re-crystallisation from solvent: iso-pentanol/diethyl ether", "diethyl ether/isoamyl alcohol"),
    ("Re-crystallisation from solvent: tert-butanol/n-hexane", "hexane/t-butanol"),
    ("Re-crystallisation from solvent: Hexane/acetong", "acetone/hexane"),
    ("Re-crystallisation from solvent: aceton and ethanol", "acetone/ethanol"),
    ("Re-crystallisation from solvent: Aqueous Solutions of Hydrochloric Acid", "hydrochloric acid"),
    ("Re-crystallisation from solvent: actone/hexane", "acetone/hexane"),
    ("Re-crystallisation from solvent: methylene chloride / n-hexan", "dichloromethane/hexane"),
    ("Re-crystallisation from solvent: acetnitrile and methanol", "acetonitrile/methanol"),
    ("Re-crystallisation from solvent: CL2CHCHCl2/i-PrOH", "1,1,2,2-tetrachloroethane/2-propanol"),
])
def test_recrystallisation_prefix_with_new_synonyms_and_typos(raw, expected):
    assert parse_value(raw) == expected
