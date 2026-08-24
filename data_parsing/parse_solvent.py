import os
import sys
import re
from typing import List, NamedTuple, Optional, Tuple
import pandas as pd
from rapidfuzz import process, fuzz

from solvent_data import (
    WHITELIST as _WHITELIST,
    SYNONYMS as _SYNONYMS,
    NON_SOLVENT_KW as _NON_SOLVENT_KW,
    METHOD_PREFIXES as _METHOD_PREFIXES,
    DESCRIPTION_KW as _DESCRIPTION_KW,
)

# ============================================================
# parse_solvent.py — cleaning the 'solvent' column
#
# Philosophy (same as parse_temperature): fail-fast on unique values.
# Iteration:
#   ⛔ unknown component → decide:
#      - genuinely missing solvent → add to WHITELIST (or SYNONYMS if a variant, in solvent_data.py)
#      - noise / method            → add to NON_SOLVENT_KW or METHOD_PREFIXES (solvent_data.py)
#   Then rerun: picks up where it stopped.
#   Pass 'restart' to reanalyze everything from scratch.
#
# Differences vs the template:
#   - Output = canonical string (not a float)
#   - Iterates over UNIQUE VALUES sorted by frequency (not row by row)
#   - Whitelist + synonyms instead of format regexes
# ============================================================

SOURCE_FILE = '../csd_all.csv'
COLUMN = 'solvent'
CLEAN_COLUMN = 'solvent_clean'

DESCRIPTION_MAX_LEN = 80   # beyond this → free-form description → None, no stop
BATCH_SIZE = 5             # number of unknown components shown per main() run


class Diagnosis(NamedTuple):
    """Diagnosis for a value resolved to `SolventParser.UNKNOWN`: why, and what to suggest."""
    preprocessed: Optional[str]
    unknown: Optional[str]
    suggestions: List[Tuple[str, float, bool]]  # (name, score, is_positional_isomer)


class SolventParser:
    """Fail-fast parser for the 'solvent' column: whitelist + synonyms.

    `parse(raw)` returns the canonical form (str), None (legitimate value
    with no solvent), or `SolventParser.UNKNOWN` (unrecognized component).
    """

    # ── QUALIFIERS TO STRIP (before splitting into components) ────────────────
    # hot/cold/dry/aqueous/XX% etc. don't change the solvent's identity.
    # Exception: the isotope-labeling group below (deuterated, perdeuterated...)
    # DOES change the identity (chloroform vs chloroform-d) but is stripped here
    # anyway — known limitation, tracked as xfail (test_deuterated_solvent_keeps_isotope_info).
    LEADING_QUALIFIER_RE = re.compile(
        r"""
        ^(
            # temperature / physical state
            very | highly | hot | heated | boiling | refluxing | cold | cooled | chilled
          | warm | dry | wet | anhydrous | anh\.\s* | absolute | abs\.
          | glacial | fuming | undistilled | dehydrated | diluted | dilute | dil\. | dil\b
          | rectified | acidic | basic | acidified | basified | chlorinated\s+ | bench\s+

            # colors
          | colorless | colourless | pale | light\s+ | dark\s+ | bright\s+ | clear\s+
          | yellow | red | green | blue | orange | purple | violet | pink | brown
          | white | black | grey | gray

            # concentration / purity
          | concentrated | concetrated | conc\. | conc\b | supersaturated | saturated
          | staturated | water-saturated | slightly\s+ | alkaline | aqueous | gaseous
          | liquid | molten | racemic-?\s* | optionally\s+ | room\s+temp\w*\s+ | readsorbed\s+

            # gas treatment
          | d\d+-
          | de-aerated | deaerated | deoxygenated | degassed | distilled | redistilled
          | freshly\ distilled

            # isotope labeling
          | perdeutero | perdeuterated
          | deuterated | ddeuterated | deutero | moist | damp | pure | fresh
          | excess\s+

            # numeric / percent prefixes
          | [><=~]?\d+[\.,]?\d*\s*mol\s*%\s*
          | [><=~]?\d+[\.,]?\d*\s*(?:[mMnN%]|percent)(?!\w)\s*
          | [><=~]?\d+[\.,]?\d*\s*%\s*
          | [><=~]?\d+\s*%\s+
          | \d+-\d+[cCkK]?\s+
          | \d+[\.,]?\d*\s+degrees?\s*
        )
        """,
        re.IGNORECASE | re.VERBOSE,
    )

    # ── UTILITY REGEXES ─────────────────────────────────────────────────────
    # Strips parenthesized annotations that aren't part of the chemical name:
    # ratios/concentrations, but also role/state/process qualifiers.
    # Rule: strip if the content starts with a digit (always a quantity),
    #       is a fraction symbol (v/v, w/w…), or matches a known qualifier.
    # Examples stripped: (1:2), (60-90), (50%), (48% in H2O), (v/v), (9:1),
    #                    (aqueous), (minor), (wet), (rt), (compound), (ligand)
    # Examples kept: (trimethylsilyl), (dimethyl), (II)
    PAREN_ANNOTATION_RE = re.compile(
        r"""
        \s*\(\s*
        (?:
            # numeric content / ratios
            \d[^)]*
          | (?:v/v|w/v|w/w|v/w|v:v|w:v|w:w|v:w)[^)]*
          | vol[^)]*
          | wt[^)]*
          | at\s+[^)]*
          | [a-z],[a-z][^)*]*

            # generic descriptors
          | aqueous
          | total\b[^)]*
          | aq\.?
          | conc\.?
          | anhydr\.?
          | anh\.?
          | [+\-±]
          | [ivxlrs]{1,4}

            # process / fraction terms
          | layer
          | antisolvent | anti-solvent
          | wash(?:ing)?
          | precipitant
          | co-?solvent
          | column[^)]*
          | eluent[^)]*
          | mobile\s+phase[^)]*
          | fraction[^)]*
          | filtrate[^)]*

            # pH / pD
          | pH\s*=\s*\d[^)]*
          | pD\s*=\s*\d[^)]*

            # qualifiers
          | minor\b | major\b | trace\b | wet\b | dry\b
          | rt\.?\b | r\.t\.?\b | room\s+temp\w*

            # range / approx
          | above\s+[^)]* | below\s+[^)]* | over\s+[^)]* | to\s+dryness
          | ca\.?\s*\d[^)]* | approx\.?\s*\d[^)]* | ~\s*\d[^)]*

            # generic chemistry nouns
          | compound[^)]* | product[^)]* | ligand[^)]* | complex[^)]*

            # fallback: short word + number
          | [a-z]{2,5}\s+\d+[^)]*
        )
        \s*\)
        """,
        re.IGNORECASE | re.VERBOSE,
    )
    LEADING_RATIO_RE = re.compile(r'^\d+\s*[:/]\s*\d+\s+')

    # ── TRAILING NOISE TO STRIP (crystallization conditions/method) ───────────
    # Removes the crystallization condition or method trailing the solvent
    # name: temperature ("at 25C", "at 238K"), method ("by slow evaporation",
    # "under nitrogen"), or a generic noun ("solution", "mixture").
    # Examples stripped: "toluene by slow cooling" -> "toluene"
    #                     "tetrahydrofuran at 238K" -> "tetrahydrofuran"
    TRAILING_NOISE_RE = re.compile(
        r"""
        \s+
        (
            under\s+\S+                                            # under nitrogen / argon / vacuum...
          | cooled\s+to(?:\s+\S+)*                                  # cooled to -20 deg
          | at\s+[\d\.\-°]+\s*\S*                                   # at 25 C / at -10 deg.
          | at\s+room\s+temp\w*                                     # at room temperature
          | by\s+slow\s+evaporation
          | by\s+evaporation
          | by\s+slow\s+cooling
          | slow\s+evaporation
          | slow\s+cooling
          | slow\s+diffusion
          | at\s+rt\.?
          | at\s+r\.t\.?
          | rt\.?\b
          | r\.t\.?\b
          | at\s+(?:low|high|ambient|elevated|reduced|room|rom|reflux)\s+(?:temp\w*|concentration\w*|pressure\w*)
          | layered\s+by\s+\S+
          | layered\s+with\s+\S+
          | as\s+reagent\s+and\s*
          | as\s+reagent
          | on\s+application\s+of\s+pressure
          | under\s+pressure
          | under\s+nitrogen
          | under\s+argon
          | under\s+inert
          | at\s+ph\s*[\d\.]+
          | at\s+(?:ca\.?\s+|approximately\s+|about\s+)?[\d\.\-°]+(?:\s*\S+)*?(?=\s*$)   # at ca. 40 (lookahead: end of string)
          | by\b
          | at\s+\d+[\.,]?\d*\s*deg
          | at\s+\d+[\.,]?\d*\s*°
          | \s*\d[\d\-\.]*\s*[kK]\b                                  # trailing temperature in K
          | \s*[-\d\.]+\s*°?[cC]\b                                   # trailing temperature in degC
          | solution\s+at\b.*
          | solution | mixture | mix | solvate | solvent
          | evaporation | diffusion | system | atmosphere
          | from\s+nmr\S*\b
          | stored\b | kept\b | layer\b | conc\.?
          | by\s+(?:freezing|melting|cooling|crystalliz\w+|crystallis\w+|sublim\w+|precipit\w+|evapor\w+)
          | \d+[\.,]?\d*\s*%(?:[-–]\d+[\.,]?\d*\s*%)?
        )
        \s*$
        """,
        re.IGNORECASE | re.VERBOSE,
    )
    TRAILING_RATIO_RE = re.compile(
        r'(\s*=\s*\d[\d\s/:.\-]*|\s+\d+(?:\.\d+)?\s*[:/]\s*\d+(?:\.\d+)?(\s*[:/]\s*\d+(?:\.\d+)?)*'
        r'|\s+\d+\s*-\s*\d+(\s*-\s*\d+)*|\s+\d+(?:\.\d+)?)\s*$'
    )

    # Sentinel — distinguishes "unknown component → stop" from a legitimate None
    UNKNOWN = object()

    def __init__(self, whitelist=_WHITELIST, synonyms=_SYNONYMS, non_solvent_kw=_NON_SOLVENT_KW,
                 method_prefixes=_METHOD_PREFIXES, description_kw=_DESCRIPTION_KW,
                 description_max_len=DESCRIPTION_MAX_LEN):
        self.whitelist = whitelist
        self.synonyms = synonyms
        self.non_solvent_kw = non_solvent_kw
        self.method_prefixes = method_prefixes
        self.description_kw = description_kw
        self.description_max_len = description_max_len

    # ── INTERNAL PIPELINE STEPS ─────────────────────────────────────────────

    def _preprocess(self, raw, trace=None):
        """Returns the cleaned (lowercase) string ready for analysis, or None.

        If `trace` is a list, each key step of the pipeline is appended to it
        (useful for debugging a specific case with `python parse_solvent.py debug "<value>"`).
        Orchestrates the named steps below; logging happens here so each
        helper can stay a plain string -> string (or -> bool) transform.
        """
        def log(label, value=...):
            if trace is not None:
                trace.append(f"{label} : {sl!r}" if value is ... else f"{label} : {value!r}")

        if pd.isna(raw):
            log('missing value (NaN)', None)
            return None
        s = str(raw).strip()
        if not s or s.lower() in ('nan', 'none', '?'):
            log('empty value/placeholder', None)
            return None
        if len(s) > self.description_max_len:
            log(f'length > {self.description_max_len} (free-form text)', None)
            return None

        sl = self._normalize_basic(s)
        sl = self._strip_leading_temperature(sl)
        log('after basic cleanup + leading temperature strip')

        # Strip a leading ratio (e.g. "1:1 molar mixture of...") before the prefix check
        sl = self.LEADING_RATIO_RE.sub('', sl).strip()

        # Remove method prefixes in a loop (e.g. "Re-cryst from: mixture of ...")
        # BEFORE the non-solvent check: a prefix like "Re-crystallisation from
        # solvent:" contains the word 'crystallisation' (present in non_solvent_kw)
        # and would otherwise be wrongly rejected before it could be stripped.
        sl, prefix_stripped = self._strip_method_prefixes(sl, max_iterations=4)
        log('after 1st pass of method prefixes' if prefix_stripped else '1st pass of method prefixes (no match)')

        sl = self._strip_leading_sample_code(sl)

        # Confirmed non-solvents + description keywords (checked after prefix strip)
        hit = self._matched_keyword(sl, self.non_solvent_kw)
        if hit is not None:
            log(f"matched non_solvent_kw {hit!r} -> None")
            return None
        hit = self._matched_keyword(sl, self.description_kw)
        if hit is not None:
            log(f"matched description_kw {hit!r} -> None")
            return None

        sl = self._strip_qualifiers_loop(sl)
        log('after stripping qualifiers (hot/dry/aqueous/95%...)')

        # 2nd pass of method prefixes: stripping qualifiers can reveal a hidden prefix
        # e.g. '2M solution HCl' → qualifier strip → 'solution HCl' → prefix strip → 'HCl'
        sl, _ = self._strip_method_prefixes(sl, max_iterations=2)
        log('after 2nd pass of method prefixes')

        sl = self._strip_ratios_and_suffixes(sl)
        log('after stripping ratios/condition suffixes (final)')

        if self._is_residual_noise(sl):
            log('residual = known noise residue -> None')
            return None
        return sl

    def _normalize_basic(self, s):
        """Lowercases and normalizes whitespace/punctuation/temperature-unit noise."""
        sl = s.lower()
        sl = sl.strip("'\"")                       # stray quotes/apostrophes at start/end
        sl = re.sub(r'[–—]', '-', sl)             # en-dash / em-dash → ASCII hyphen
        sl = re.sub(r'[\n\r\t_]+', ' ', sl)       # newlines/tabs/underscores → space
        sl = re.sub(r'  +', ' ', sl).strip()      # multiple spaces → single space

        # Normalize "Xat N" → "X at N" (concatenations with no space before 'at')
        sl = re.sub(r'([a-z])(at\s+[\d\-])', r'\1 \2', sl)
        # Normalize "30 deg.C" / "30 deg C" / "-35deg.C" → "30°C"
        sl = re.sub(r'(\d)\s*deg\.?\s*([cCkK])\b', r'\1°\2', sl)  # digit glued to it: -35deg.C
        sl = re.sub(r'\bdeg\.?\s*([cCkK])\b', r'°\1', sl)          # with space or at word start
        # Normalize "N, N-dimethyl" → "N,N-dimethyl" (avoids a spurious split on comma+space)
        sl = re.sub(r'\b([nNoO]),\s+([nNoO])-', r'\1,\2-', sl)
        # Normalize spaces around / (e.g. "v/ v" → "v/v", "acetone/ water" → "acetone/water")
        sl = re.sub(r'\s*/\s*', '/', sl)
        # Normalize "n- hexane" → "n-hexane" (space after the 'n-' abbreviation's hyphen)
        sl = re.sub(r'\bn-\s+', 'n-', sl)
        return sl

    def _strip_leading_temperature(self, s):
        """Strips a leading bare temperature (e.g. '-35 c from ether', '120 celsius activated')."""
        return re.sub(r'^-?\d+\.?\d*\s*(?:°?[ckCK]\b|celsius\b|kelvin\b|degrees?\b)\s*', '', s, flags=re.IGNORECASE).strip()

    def _strip_method_prefixes(self, s, max_iterations):
        """Repeatedly strips a leading ratio then a recognized method prefix
        (e.g. "Re-cryst from: mixture of ..."). Returns (new_s, stripped_any).
        """
        stripped_any = False
        for _ in range(max_iterations):
            s = self.LEADING_RATIO_RE.sub('', s).strip()  # re-strip after each prefix removed
            stripped_this = False
            for prefix in self.method_prefixes:
                if s.startswith(prefix):
                    remainder = s[len(prefix):]
                    # Avoid cutting mid-word (e.g. "from a" latching onto "from acetone")
                    if prefix and prefix[-1].isalpha() and remainder and remainder[0].isalpha():
                        continue
                    s = remainder.strip(' :,')
                    stripped_any = True
                    stripped_this = True
                    break
            if not stripped_this:
                break
        return s, stripped_any

    def _strip_leading_sample_code(self, s):
        """Strips a leading sample code (e.g. 'mds088 from a mixture of...') before from/of."""
        return re.sub(r'^[a-z]{1,5}\d{1,4}\s+(?=(?:from|of)\b)', '', s).strip()

    @staticmethod
    def _matched_keyword(s, keywords):
        """Returns the first keyword from `keywords` found as a substring of `s`, or None."""
        return next((kw for kw in keywords if kw in s), None)

    def _strip_qualifiers_loop(self, s):
        """Strips leading qualifiers (hot/dry/aqueous/95%...) in a loop, also
        trimming a residual leading hyphen/period after each pass. Differs from
        `_strip_qualifiers` (used by `_resolve`) by that extra lstrip and by
        bailing out if a pass would empty the string.
        """
        for _ in range(3):
            new = self.LEADING_QUALIFIER_RE.sub('', s).strip().lstrip('-. ')
            if new == s or not new:
                break
            s = new
        return s

    def _strip_ratios_and_suffixes(self, s):
        """Strips parenthesized ratios, condition suffixes (at -30°C, by slow
        evaporation...), and trailing ratio notations (1:1, v/v, 40/60...).
        """
        s = self.PAREN_ANNOTATION_RE.sub('', s).strip()
        s = self.LEADING_RATIO_RE.sub('', s).strip()
        # Strip ratio annotations = N/M or = N:M in the middle/end of the string
        s = re.sub(r'\s*=\s*[\d]+\s*[/:][\d]+(\s*[/:][\d]+)*', '', s).strip()
        # Strip in a loop: several suffixes can chain together
        for _ in range(5):
            new = self.TRAILING_NOISE_RE.sub('', s).strip().rstrip('.').strip()
            if new == s:
                break
            s = new

        # Strip N/M ratios at the end of the string (e.g. "ethanol 1/1", "petroleum ether 40/60")
        s = re.sub(r'\s+\d+/\d+(/\d+)*\s*$', '', s).strip()
        # Strip N:M (v/v) ratios outside parentheses at the end of the string (with or without a leading hyphen/comma)
        s = re.sub(r'\s*[-–]\s*\d+\s*:\s*\d+(\s*:\s*\d+)*\s*(?:v/v|w/w|v/w|w/v)?\s*$', '', s).strip()
        s = re.sub(r'\s+\d+\s*:\s*\d+(\s*:\s*\d+)*\s*(?:v/v|w/w|v/w|w/v)?\s*$', '', s).strip()
        s = re.sub(r',\s*\d+\s*[:/]\s*\d+(\s*[:/]\s*\d+)*\s*(?:v/v|w/w|v/w|w/v)?\s*$', '', s).strip()
        s = re.sub(r'\s+(?:v/v|w/w|v/w|w/v)\s*$', '', s).strip()
        s = re.sub(r'\s+\d+/\d+(/\d+)*\s*$', '', s).strip()   # 2nd pass after stripping v/v
        # Residual trailing hyphen/comma (after stripping the ratio)
        s = re.sub(r'\s*[-–]\s*$', '', s).strip()
        s = s.rstrip(',').strip()
        return s

    def _is_residual_noise(self, s):
        """True if the leftover string is known noise with no solvent information."""
        return bool(
            not s or s in (
                'the', 'a', 'an', 'of', 'at', 'by', 'in', 'on', 'to', 'from', 'with',
                'solvent', 'mixture', 'oil', '-', '?', 'no', 'v', 'acid', 'base', 'salt',
                'solution', 'solutions', 'heated', 'cooled', 'stirred', 'filtered',
                'freezing', 'melting', 'boiling', 'sublimation', 'degree', 'degrees',
            ) or s.isdigit() or re.fullmatch(r'[?!\-.*_\s]+', s)
        )

    def _strip_qualifiers(self, s):
        """Strips leading qualifiers (hot, dry, aqueous, 95%…) in a loop."""
        for _ in range(3):
            new = self.LEADING_QUALIFIER_RE.sub('', s).strip()
            if new == s:
                break
            s = new
        return s

    def _split_components(self, s):
        """Splits a solvent string into individual components."""
        # Hyphen separator before N,N- / N,O- etc. (IUPAC substituent names)
        s = re.sub(r'(?<=[a-z])-([no],[no]-)', r'/\1', s, flags=re.IGNORECASE)
        s = re.sub(r'\s*\\\s*', '/', s)        # backslash → /
        s = re.sub(r'\s*/\s*', '/', s)        # spaces around /
        s = re.sub(r'\s+and\s*/', '/', s)     # 'and/' (and stuck to the slash)
        s = re.sub(r'\s*\blayered?\s+(?:over|with|on|onto|by)?\s*', '/', s)  # 'layered over/with/on/onto/bare'
        s = re.sub(r'(?<=[a-z]) layer (?=[a-z])', '/', s, flags=re.IGNORECASE)  # 'thf layer hexane' → 'thf/hexane'
        s = re.sub(r'\s+-\s+', '/', s)        # hyphen with spaces on both sides
        s = re.sub(r'\s+a?and\s+', '/', s)    # 'and' (typo 'aand' included)
        s = re.sub(r'\s*&\s*', '/', s)        # '&'
        s = re.sub(r'\s*\+\s*', '/', s)      # '+'
        s = re.sub(r'\s+or\s+', '/', s)       # 'or'
        s = re.sub(r'\s+containing\s+', '/', s)  # 'acetone containing water' → two components
        s = re.sub(r'\s+followed\s+by\s+', '/', s)  # 'ethanol followed by acetone' → two components
        s = re.sub(r'\s+on(?:to)?\s+', '/', s)   # 'Et2O on/onto MeOH' (layered) → two components
        s = re.sub(r'\s+over\s+', '/', s)   # 'hexane over THF' (layering) → two components
        s = re.sub(r'\s+to\s+', '/', s)     # 'hexanes to methylene chloride' → mixture separator
        s = re.sub(r'\s+vs\.?\s+', '/', s)  # 'acetone vs hexane' → separator
        s = re.sub(r'\s*;\s*', '/', s)       # semicolon separator
        s = re.sub(r'(?<=[a-z])\s*:\s*(?=[a-z\d])', '/', s)  # colon separator between letters (or letter:digit)
        s = re.sub(r'(?<=[a-z]{5})\.\s+(?=[a-z])', '/', s)  # period after 5+ letters = separator
        s = re.sub(r',\s+', '/', s)           # comma+space (not in IUPAC names)
        s = re.sub(r'(?<=[a-z]{2}),(?=[a-z]{2})', '/', s)  # comma with no space between 2 words (e.g. dcm,meoh)

        if '/' in s:
            return [c.strip() for c in s.split('/') if c.strip()]

        # Hyphen as separator: 2+ letters before AND 2+ letters after (avoids -d6, -d4…)
        # Only if every part resolves to a known solvent
        # (avoids breaking "di-isopropyl ether", "n-butyl alcohol"…)
        parts = re.split(r'(?<=[a-z][a-z])-(?=[a-z][a-z])', s)
        if len(parts) > 1:
            if all(self._resolve(p.strip()) in self.whitelist for p in parts if p.strip()):
                return [c.strip() for c in parts if c.strip()]

        # Fallback: try every hyphen position (e.g. "chloroform-n-hexane")
        if '-' in s:
            positions = [m.start() for m in re.finditer(r'-', s)]
            for pos in positions:
                left, right = s[:pos], s[pos+1:]
                if not left or not right:
                    continue
                rl = self._resolve(left)
                rr = self._resolve(right)
                if rl in self.whitelist and rr in self.whitelist:
                    return [left, right]

        # Last resort: space as separator — only if every chunk
        # resolves to a known solvent (avoids breaking "diethyl ether", etc.)
        words = s.split()
        if len(words) >= 2:
            # 2 parts
            for i in range(1, len(words)):
                left  = ' '.join(words[:i])
                right = ' '.join(words[i:])
                if self._resolve(left) in self.whitelist and self._resolve(right) in self.whitelist:
                    return [left, right]
            # 3 parts
            for i in range(1, len(words) - 1):
                for j in range(i + 1, len(words)):
                    p1 = ' '.join(words[:i])
                    p2 = ' '.join(words[i:j])
                    p3 = ' '.join(words[j:])
                    if all(self._resolve(p) in self.whitelist for p in [p1, p2, p3]):
                        return [p1, p2, p3]

        return [s.strip()]

    def _resolve(self, component):
        """Applies qualifiers + synonym lookup, returns the canonical form or the original string."""
        # Direct synonym lookup before any stripping (avoids the qualifier-strip wiping out the component)
        if component in self.synonyms:
            return self.synonyms[component]
        # A bare temperature (e.g. '-20°c', '343k') → a condition, not a solvent
        if re.fullmatch(r'-?\d+\.?\d*\s*°?\s*[ckCK]', component):
            return None
        # A bare number or ratio (e.g. '1.6', '12', '100:15', '1:2:1') → ignored
        if re.fullmatch(r'[\d]+(?:[.,]\d+)?(?:\s*[:/]\s*[\d]+(?:[.,]\d+)?)*', component):
            return None

        c = self._strip_resolve_ratios_and_parens(component)
        c = self._strip_resolve_trailing_keywords(c)

        if self._is_resolve_noise(c):
            return None

        # Strip 'from'/'of'/'containing' / 'N drops of' / 'a little' prefixes at the start of an individual component
        c2 = self._strip_resolve_prefix(c)
        if c2 != c:
            c2 = self._strip_qualifiers(c2).lstrip('-. ')
            if not c2:
                return None
            result = self.synonyms.get(c2, c2)
            if result in self.whitelist or result is None:
                return result

        resolved = self._unwrap_own_abbreviation(c)
        if resolved is not None:
            return resolved
        return self.synonyms.get(c, c)

    def _strip_resolve_ratios_and_parens(self, c):
        """Strips qualifiers, leading/trailing/parenthesized ratio annotations,
        and unwraps a component entirely wrapped in parentheses.
        """
        c = self._strip_qualifiers(c)
        c = c.lstrip('-. ')   # residual hyphen/period after stripping qualifiers (e.g. perdeutero-X)
        c = self.LEADING_RATIO_RE.sub('', c).strip()     # strips "12:1 " at the start of a component (e.g. after splitting on /)
        c = self.TRAILING_RATIO_RE.sub('', c).strip()   # strips "1:1", "1:2:1" at the end of a component
        c = self.PAREN_ANNOTATION_RE.sub('', c).strip()  # strips "(1:1)", "(v/v)" etc.
        if c.count(')') > c.count('('):              # unmatched trailing ')' (e.g. 'ethyl acetate)')
            c = c.rstrip(')')
        # Component entirely wrapped in parentheses (e.g. "solvent (acetonitrile)" → "(acetonitrile)")
        # → unwrap, the content is the real solvent name
        m = re.fullmatch(r'\((.*)\)', c)
        if m and '(' not in m.group(1) and ')' not in m.group(1):
            c = m.group(1).strip()
        return c

    def _strip_resolve_trailing_keywords(self, c):
        """Strips trailing method/condition keywords left on a component
        (solution, mixture, binary, co-solvent, by..., was..., as..., a
        residual isolated letter or comma, leftover temperature/duration,
        saturated, ionic liquid, containing, reaction).
        """
        # Strip trailing keywords on the component (applied before solution/mixture, for compound strings)
        c = re.sub(r'\s+\bto\b\s*$', '', c).strip()             # 'diethyl ether solution to' → 'diethyl ether solution'
        c = re.sub(r'\s+solutions?\s*$', '', c).strip()  # '1,2-difluorobenzene solution(s)' → '1,2-difluorobenzene'
        c = re.sub(r'[-]mixture\s*$', '', c, flags=re.IGNORECASE).strip()  # 'hexane-mixture' → 'hexane'
        c = re.sub(r'\s+binary\b.*$', '', c, flags=re.IGNORECASE).strip()  # 'water binary' → 'water'
        c = re.sub(r'\s+co-?solvents?(?:\s+system)?\s*\.?\s*$', '', c, flags=re.IGNORECASE).strip()  # 'ethanol co-solvent system.' → 'ethanol'
        c = re.sub(r'\s+by\w*(?:\s+\S+)*\s*$', '', c, flags=re.IGNORECASE).strip()  # 'hexane bydiffusion' / 'dcm by --' → solvent
        c = re.sub(r'\s+was\s+\w+\s*$', '', c, flags=re.IGNORECASE).strip()  # 'naphthalene was present' → 'naphthalene'
        c = re.sub(r'\s+as\s+\w[\w\s]*$', '', c, flags=re.IGNORECASE).strip()  # 'antimony pentachloride as reagent' → 'antimony pentachloride'
        c = re.sub(r'\s+[a-z]\s*$', '', c).strip()    # residual isolated letter (e.g. 'dichloromethane a')
        c = c.rstrip(',').strip()                       # trailing comma (e.g. 'methanol,')
        # Strip remaining temperature/duration conditions from the component after splitting on /
        c = re.sub(r'\s+for\s+\d+.*$', '', c).strip()           # 'pentane for 8 days' → 'pentane'
        c = re.sub(r'\s+[-−]?\d+\.?\d*\s*°?[cCkK]\b.*$', '', c).strip()  # 'pentane -30°c ...' → 'pentane'
        c = re.sub(r'\s+mixture\b.*$', '', c).strip()            # 'hexane mixture via slow ...' → 'hexane'
        c = re.sub(r'\s+mixed\b.*$', '', c).strip()              # 'heptane mixed 1:1' → 'heptane'
        c = re.sub(r'\s+at\s+(?:rt|r\.t\.?|room\s+temp\w*).*$', '', c, flags=re.IGNORECASE).strip()
        c = re.sub(r'\s+\bsaturated\b.*$', '', c, flags=re.IGNORECASE).strip()  # 'benzene saturated' → 'benzene'
        c = re.sub(r'\s+ionic\s+liquids?\b.*$', '', c, flags=re.IGNORECASE).strip()  # 'bmim acetate ionic liquid' → 'bmim acetate'
        c = re.sub(r'\s+contain(?:ing|g)\b.*$', '', c, flags=re.IGNORECASE).strip()  # 'etoh containg nanocluster 3' → 'etoh' (typo 'containg' included)
        c = re.sub(r'\s+reaction\b.*$', '', c, flags=re.IGNORECASE).strip()  # 'methanol reaction' → 'methanol'
        return c

    @staticmethod
    def _is_resolve_noise(c):
        """True if `c` is leftover noise (empty, a residual bare temperature, an
        isolated preposition, or an imidazolium/pyrrolidinium ionic liquid)
        rather than an actual solvent component.
        """
        if not c:
            return True
        # A bare temperature after stripping (e.g. '85°c' left over from 'heated to 85°c for 20 min')
        if re.fullmatch(r'-?\d+\.?\d*\s*°?\s*[ckCK]', c):
            return True
        # Residue after stripping a qualifier or splitting on ' and ': an isolated preposition/condition
        if re.match(r'^(?:to|until|while|before|after|upon|via|using|with|without|over|under|through)\b', c, re.IGNORECASE):
            return True
        # Imidazolium salts / ionic liquids → not in the whitelist
        if 'imidazolium' in c or 'pyrrolidinium' in c:
            return True
        return False

    @staticmethod
    def _strip_resolve_prefix(c):
        """Strips 'from'/'of'/'containing' / 'N drops of' / 'a little' prefixes
        at the start of an individual component.
        """
        return re.sub(r'^(?:a\s+little\s+|little\s+|a\s+small\s+amount\s+of\s+|an?\s+(?=\w)|from|of|containing|\d+\s+drops?\s+of\s*|a?\s*(?:few|one|two|several|some)?\s*drops?\s+of)\s*', '', c).strip()

    def _unwrap_own_abbreviation(self, c):
        """If `c` is "full name (own abbreviation)", returns the resolved
        canonical name when that resolution is in the whitelist (the
        abbreviation is just a reminder, e.g. "n,n-diethylformamide (def)");
        otherwise returns None.
        """
        m = re.match(r'^(.*\S)\s*\([a-z0-9]{2,8}\)$', c)
        if not m:
            return None
        candidate = m.group(1).strip()
        resolved = self.synonyms.get(candidate, candidate)
        return resolved if resolved in self.whitelist else None

    # ── PUBLIC API ───────────────────────────────────────────────────────────

    def parse(self, raw):
        """
        Returns the canonical form (str), None (legitimate value with no solvent),
        or `SolventParser.UNKNOWN` (unrecognized component → the caller should stop).

        If `raw` is a list, returns a list of results (one per element).
        """
        if isinstance(raw, list):
            return [self.parse(v) for v in raw]

        sl = self._preprocess(raw)
        if sl is None:
            return None

        # Direct whitelist match (before splitting — avoids breaking hyphenated names)
        if sl in self.whitelist:
            return sl

        # Global synonym lookup on the whole string
        resolved = self.synonyms.get(sl)
        if resolved is not None:
            return resolved

        components = self._split_components(sl)
        if not components:
            return None

        normalized = []
        for c in components:
            canonical = self._resolve(c)
            if canonical is None:
                continue  # non-solvent additive (e.g. electrolyte salt) → ignored
            if canonical not in self.whitelist:
                # Re-split the component (hyphen or space) — case after splitting on /
                sub = self._split_components(c)
                if len(sub) > 1:
                    sub_ok = True
                    sub_names = []
                    for sc in sub:
                        rc = self._resolve(sc)
                        if rc is None:
                            continue
                        if rc not in self.whitelist:
                            sub_ok = False
                            break
                        sub_names.append(rc)
                    if sub_ok and sub_names:
                        normalized.extend(sub_names)
                        continue
                return self.UNKNOWN
            normalized.append(canonical)

        if not normalized:
            return None
        return '/'.join(sorted(set(normalized)))

    def diagnose(self, raw):
        """Full diagnosis for a value resolved to `UNKNOWN`: string after
        preprocessing, the offending component, and whitelist suggestions
        (with a positional-isomer warning). Meant for external callers
        (e.g. `main()`) that shouldn't have to touch the pipeline's private methods.
        """
        preprocessed = self._preprocess(raw)
        unknown = self._find_unknown_component(raw)
        suggestions = [
            (name, score, self._is_positional_isomer(unknown or '', name))
            for name, score in self._fuzzy_suggest(unknown or '')
        ]
        return Diagnosis(preprocessed, unknown, suggestions)

    def _find_unknown_component(self, raw):
        """Returns the first unrecognized component (for the stop message)."""
        sl = self._preprocess(raw)
        if sl is None:
            return None
        if sl in self.synonyms:
            return None
        for c in self._split_components(sl):
            r = self._resolve(c)
            if r is not None and r not in self.whitelist:
                return c
        return None

    @staticmethod
    def _is_positional_isomer(a, b):
        """True if a and b differ only by a position digit (e.g. 1,2- vs 1,3-)."""
        a_norm = re.sub(r'\d', 'N', a)
        b_norm = re.sub(r'\d', 'N', b)
        return a_norm == b_norm and a != b

    def _fuzzy_suggest(self, component, n=5):
        """Top n whitelist candidates by similarity (fuzz.ratio)."""
        hits = process.extract(component, list(self.whitelist), scorer=fuzz.ratio, limit=n)
        return [(name, score) for name, score, _ in hits]

    def debug(self, raw):
        """Prints the parsing pipeline step by step for a given value.

        Usage: python parse_solvent.py debug "<value>"
        """
        print(f"RAW          : {raw!r}")

        trace = []
        sl = self._preprocess(raw, trace=trace)
        print("\n-- _preprocess --")
        for step in trace:
            print(f"  {step}")
        print(f"PREPROCESSED : {sl!r}")

        if sl is None:
            print("\nRESULT       : None  (legitimate, not a solvent)")
            return

        if sl in self.whitelist:
            print("\n-- direct whitelist match (before split) --")
            print(f"RESULT       : {sl!r}")
            return

        resolved = self.synonyms.get(sl)
        if resolved is not None:
            print(f"\n-- global synonym on the whole string : {sl!r} -> {resolved!r} --")
            print(f"RESULT       : {resolved!r}")
            return

        components = self._split_components(sl)
        print(f"\n-- _split_components --\n  {components}")

        print("\n-- _resolve per component --")
        normalized = []
        for c in components:
            canonical = self._resolve(c)
            status = 'OK' if canonical in self.whitelist else ('ignored (non-solvent)' if canonical is None else 'UNKNOWN')
            print(f"  resolve({c!r}) -> {canonical!r}  [{status}]")
            if canonical is None:
                continue
            if canonical in self.whitelist:
                normalized.append(canonical)

        result = self.parse(raw)
        print(f"\nRESULT       : {'UNKNOWN (stop fail-fast)' if result is self.UNKNOWN else result!r}")
        if result is self.UNKNOWN:
            diag = self.diagnose(raw)
            print(f"  unknown component : {diag.unknown!r}")
            print("  suggestions       :")
            for name, score, is_isomer in diag.suggestions:
                warning = "  ⚠ positional isomer — check this!" if is_isomer else ""
                print(f"    {score:3.0f}  {name}{warning}")


# ── Default instance + module-level aliases (CLI / test compatibility) ───────

_default_parser = SolventParser()
_UNKNOWN = SolventParser.UNKNOWN


def parse_value(raw):
        return _default_parser.parse(raw)


def debug_parse(raw):
    return _default_parser.debug(raw)


def _progress_file():
    return f'.progress_{COLUMN}.txt'


# ── MAIN ──────────────────────────────────────────────────────────────────────

def main():
    parser = _default_parser
    restart = 'restart' in sys.argv[1:]
    progress = _progress_file()

    start = 0
    if not restart and os.path.exists(progress):
        try:
            start = int(open(progress).read().strip())
        except ValueError:
            start = 0

    print(f"Reading column '{COLUMN}'...")
    col = pd.read_csv(SOURCE_FILE, usecols=[COLUMN], low_memory=False)[COLUMN]

    # Work on unique values sorted by descending frequency:
    # the most impactful unknowns surface first.
    vc = col.dropna().value_counts()
    unique_vals = list(vc.index)
    n = len(unique_vals)

    if start >= n:
        start = 0

    print(f"Analyzing {n:,} unique values starting at index {start}...")

    errors = []   # list of (i, raw, diag) — diag: Diagnosis

    for i in range(start, n):
        raw = unique_vals[i]

        result = parser.parse(raw)
        if result is None:
            continue

        if result is SolventParser.UNKNOWN:
            diag = parser.diagnose(raw)
            # Deduplicate: no point showing the same unknown component twice
            if diag.unknown and any(e[2].unknown == diag.unknown for e in errors):
                continue
            errors.append((i, raw, diag))
            if len(errors) >= BATCH_SIZE:
                break

    if errors:
        # Checkpoint on the first unknown (resumes here on next run)
        with open(progress, 'w') as f:
            f.write(str(errors[0][0]))

        for idx, (i, raw, diag) in enumerate(errors, 1):
            print(f"\n⛔  Unknown component {idx}/{len(errors)}  (value #{i + 1}/{n},  {vc[raw]:,} occurrences)")
            print(f"   Raw          : {raw!r}")
            if diag.preprocessed and diag.preprocessed != raw.lower().strip():
                print(f"   After strip  : {diag.preprocessed!r}")
            print(f"   Unknown      : {diag.unknown!r}")
            print(f"\n   Whitelist suggestions (fuzz.ratio) :")
            for name, score, is_isomer in diag.suggestions:
                warning = "  ⚠ positional isomer — check this!" if is_isomer else ""
                print(f"     {score:3.0f}  {name}{warning}")
        print()
        print("→ Solvent missing from the whitelist  : add to WHITELIST (or SYNONYMS if a variant, in solvent_data.py)")
        print("→ Spelling variant of an existing one : add to SYNONYMS")
        print("→ Noise / method                      : add to NON_SOLVENT_KW or DESCRIPTION_KW")
        print("→ Then rerun (resumes here). Pass 'restart' to reanalyze everything.")
        sys.exit(1)

    # ── Everything recognized → save ──────────────────────────────────────────
    if os.path.exists(progress):
        os.remove(progress)

    print("✅ All components recognized. Reading the full file and saving...")
    df = pd.read_csv(SOURCE_FILE, low_memory=False)

    mapping = {}
    for v in vc.index:
        r = parser.parse(v)
        mapping[v] = None if (r is SolventParser.UNKNOWN or r is None) else r

    results = df[COLUMN].map(lambda x: mapping.get(x) if pd.notna(x) else None)

    if CLEAN_COLUMN in df.columns:
        df = df.drop(columns=[CLEAN_COLUMN])
    idx = df.columns.get_loc(COLUMN) + 1
    df.insert(idx, CLEAN_COLUMN, results)
    df.to_csv(SOURCE_FILE, index=False)
    print(f"✅ Column '{CLEAN_COLUMN}' added to {SOURCE_FILE}.")


if __name__ == '__main__':
    if len(sys.argv) >= 3 and sys.argv[1] == 'debug':
        debug_parse(sys.argv[2])
    else:
        main()
