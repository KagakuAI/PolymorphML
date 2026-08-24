# Solvent Parser — Maintainer Guide

`parse_solvent.py` turns the free-text `solvent` column of `csd_all.csv` into a
clean `solvent_clean` column. It works by **whitelist + synonyms**, not by trying
to parse arbitrary chemistry text — so it deliberately stops the moment it meets
a solvent it doesn't recognize, instead of silently guessing. This guide is for
whoever updates it later, most likely because new CSD entries introduced solvents
the script has never seen.

## Quick start

```bash
cd data_parsing
python parse_solvent.py                  # resumes from where it last stopped
python parse_solvent.py restart          # re-analyzes every unique value from scratch
python parse_solvent.py debug "<value>"  # step-by-step trace for one value
```

Without `restart`, the script resumes from `.progress_solvent.txt` (the index of
the first unresolved value, in a list sorted by descending frequency — so the
most impactful unknowns surface first). Delete that file, or pass `restart`, to
start over.

## What you'll see when it stops

```
⛔  Unknown component 1/5  (value #842/18219,  37 occurrences)
   Raw          : 'recrystallized from xyz-ether'
   After strip  : 'xyz-ether'
   Unknown      : 'xyz-ether'

   Whitelist suggestions (fuzz.ratio) :
     78  diethyl ether
     ...
```

For each unknown component, decide which of these applies, fix `solvent_data.py`
accordingly, then rerun — the script resumes exactly where it stopped.

| Situation | Fix |
|---|---|
| Genuinely new solvent, not yet in the whitelist | Add it to `WHITELIST` |
| Spelling/abbreviation/typo of a solvent already in `WHITELIST` | Add it to the right group in `SYNONYM_GROUPS` |
| Additive that isn't a solvent (salt, catalyst, ionic liquid…) | Add it to `NON_SOLVENT_SYNONYMS` (resolves to `None`, doesn't stop) |
| Crystallization method/noise, not a chemical name at all | Add to `NON_SOLVENT_KW` or `DESCRIPTION_KW` |
| A recurring boilerplate prefix ("Recrystallized from a mixture of...") | Add to `METHOD_PREFIXES` |

## Where things live in `solvent_data.py`

- **`WHITELIST`** — the set of canonical names. Only names in this set can appear
  in `solvent_clean`. Keep them lowercase, and prefer the common/short form over
  the IUPAC name when there's a choice already in use (e.g. `dmpu`, not the
  spelled-out heterocycle name) — check existing entries for the established style.
- **`SYNONYM_GROUPS`** — `{canonical: [variant1, variant2, ...]}`. Add the new
  spelling/abbreviation to the list under its canonical solvent.
  **Never edit `SYNONYMS` directly** — it's generated from `SYNONYM_GROUPS` +
  `NON_SOLVENT_SYNONYMS` at import time, and raises `ValueError` if a variant ends
  up assigned to two different canonicals (this is the safety net against a
  silent-overwrite bug found in an earlier version of this file).
- **`NON_SOLVENT_SYNONYMS`** — flat list of variant names that should resolve to
  `None` (additives, salts, ionic liquids — anything that isn't a solvent but
  isn't pure noise either).
- **`NON_SOLVENT_KW`** — substrings checked against the whole preprocessed string;
  if found, the value resolves to `None` without stopping. For crystallization
  descriptions that contain no extractable solvent name at all.
- **`METHOD_PREFIXES`** — boilerplate prefixes stripped before analysis
  (`"recrystallised from solvent:"`, etc.), kept sorted longest-first so a long
  prefix can't be shadowed by a shorter one that's also a substring of it.
- **`DESCRIPTION_KW`** — like `NON_SOLVENT_KW`, but checked later in the pipeline
  (after prefixes are stripped) for free-text description markers.

## Things to watch for

- **Positional isomers.** `debug`'s suggestion list flags
  `⚠ positional isomer — check this!` when a suggestion differs from the unknown
  value only by a position number (e.g. `1,2-X` vs `1,3-X`). These are different
  molecules — don't add a synonym across isomers even if the fuzzy score is high
  (80-90%+).
- **Don't strip chemical parentheses blindly.** `(trimethylsilyl)`, `(dimethyl)`,
  `(II)` are part of the name; only ratio/role/state annotations get stripped —
  see the comment above `PAREN_ANNOTATION_RE` in `parse_solvent.py` for the exact rule.
- **Single occurrence, no isomer specified** (e.g. plain `"nitropropane"`): map it
  to the most common isomer (`1-nitropropane`) rather than leaving it unresolved,
  following the existing convention — check occurrence count and look for similar
  past decisions in `WHITELIST`/`SYNONYM_GROUPS` before doing this.
- **Use `python parse_solvent.py debug "<value>"`**, not ad-hoc `python -c "..."`
  snippets, to inspect why a value resolves a certain way — it traces every
  pipeline step (`_preprocess`, `_split_components`, `_resolve` per component) and
  prints fuzzy suggestions on failure.

## Regression tests

`test_parse_solvent.py` covers known edge cases, including
`@pytest.mark.xfail` tests that document **known unfixed bugs** (e.g.
`"aqueous ethanol"` loses the water component; deuterated-solvent labels aren't
normalized to their `-d6`/`-d4`/... suffix). Run it before and after any change:

```bash
python -m pytest test_parse_solvent.py -q
```

If you fix one of the `xfail` bugs, remove its `xfail` marker — otherwise a
regression on that exact case would go unnoticed.

## After any change

1. Run the failing case through `debug` to confirm the fix.
2. `python -m pytest test_parse_solvent.py -q` — no new failures (the `xfail`
   count may drop if you fixed a known bug, but nothing should newly fail).
3. `python parse_solvent.py restart` — should finish with
   `✅ All components recognized.` and the saved column should cover every row.
