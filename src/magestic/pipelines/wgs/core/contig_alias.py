"""Bijective contig-name re-representation for TRACE's positional variant-id parse (GATE 3).

WHY THIS IS AN ADAPTER AND NOT AN ALGORITHM CHANGE
--------------------------------------------------
`trace_crispr.wgs.donor_outcome.precall_sample` builds its own variant id as
``f"{chrom}_{pos}_{ref}_to_{alt}"`` and then re-parses it with ``split("_")``,
indexing **positionally**: ``sp[1]`` is the position, ``sp[2]`` the ref and
``sp[4]`` the alt. Any underscore inside ``chrom`` shifts every index, so
``int(sp[1])`` raises ``ValueError`` mid-loop, after a partial output is on disk.

Every reference the WGS reference-mapping policy mandates carries underscore-bearing
contigs — ``blank_MAGESTIC3_pF835_guide_donor_plasmid``, ``native_2micron``,
``unplaced_contig_34``, ``plasmid_{oligo}``, ``integrated_plasmid_{oligo}`` — so the
stage is incompatible with the reference standard rather than with one odd dataset.

The permanent fix still belongs in TRACE: it should carry structured fields instead of
re-parsing a string it just formatted (a right-anchored parse would also do, since
``ref``/``alt`` are nucleotide strings and ``to`` is a literal, so
``alt=sp[-1], ref=sp[-3], pos=sp[-4], chrom="_".join(sp[:-4])`` is unambiguous).
TRACE is pinned by sha in the ``magestic[wgs]`` extra, so this module does **not**
change TRACE's behaviour. It changes only the *representation* magestic hands across
the boundary, which is the same class of adapter work ``core/vcf_adapter.py`` already
owns (remapping AD from FORMAT slot 4 to slot 2). TRACE sees contig names with no
underscore, so its positional parse is correct on its own terms, and magestic maps the
names back on the way out.

GUARANTEES
----------
* **Bijective.** ``decode(encode(name)) == name`` for every name supplied, asserted at
  construction. A collision raises rather than silently merging two contigs.
* **Underscore-free aliases.** Aliases match ``^[A-Za-z0-9]+$``, so they survive any
  positional ``split("_")``.
* **Order-independent and deterministic.** Aliases are assigned over the *sorted* input,
  so the same reference always yields the same mapping and two runs are comparable.
* **Identity for safe names.** A contig without an underscore is left untouched, so
  existing outputs on underscore-free references are byte-identical.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping, MutableMapping
from dataclasses import dataclass, field

__all__ = ["ContigAlias", "ALIAS_PREFIX", "alias_rows"]

ALIAS_PREFIX = "ctg"
_SAFE = re.compile(r"^[A-Za-z0-9]+$")


@dataclass(frozen=True)
class ContigAlias:
    """A bijective map from reference contig names to underscore-free aliases.

    Attributes
    ----------
    encode_map : dict
        Original name -> alias. Contains an entry for every supplied name, including
        the ones that needed no change (identity), so callers never have to decide.
    decode_map : dict
        Alias -> original name.
    aliased : tuple of str
        Just the names that actually required an alias, sorted. Empty tuple means the
        reference was already safe and this adapter is a no-op.
    """

    encode_map: Mapping[str, str] = field(default_factory=dict)
    decode_map: Mapping[str, str] = field(default_factory=dict)
    aliased: tuple[str, ...] = ()

    @classmethod
    def from_names(cls, names: Iterable[str]) -> ContigAlias:
        unique = sorted({str(n) for n in names if str(n) != ""})
        enc: dict[str, str] = {}
        dec: dict[str, str] = {}
        needed: list[str] = []
        # Reserve the safe names first so an alias can never shadow a real contig.
        reserved = {n for n in unique if _SAFE.match(n)}
        counter = 0
        for name in unique:
            if _SAFE.match(name):
                alias = name
            else:
                while True:
                    counter += 1
                    alias = f"{ALIAS_PREFIX}{counter:04d}"
                    if alias not in reserved and alias not in dec:
                        break
                needed.append(name)
            if alias in dec and dec[alias] != name:
                raise ValueError(
                    f"contig alias collision: {alias!r} claimed by both "
                    f"{dec[alias]!r} and {name!r}"
                )
            enc[name] = alias
            dec[alias] = name
        obj = cls(encode_map=enc, decode_map=dec, aliased=tuple(sorted(needed)))
        obj._assert_bijective()
        return obj

    def _assert_bijective(self) -> None:
        if len(self.encode_map) != len(self.decode_map):
            raise ValueError(
                f"contig alias is not bijective: {len(self.encode_map)} names map to "
                f"{len(self.decode_map)} aliases"
            )
        for name, alias in self.encode_map.items():
            if self.decode_map[alias] != name:
                raise ValueError(f"contig alias does not round-trip for {name!r}")
            if "_" in alias:
                raise ValueError(
                    f"alias {alias!r} for {name!r} still contains '_', which is the "
                    "whole point of this adapter"
                )

    def encode(self, name: str) -> str:
        """Original contig name -> underscore-free alias (identity if already safe)."""
        return self.encode_map.get(str(name), str(name))

    def decode(self, alias: str) -> str:
        """Alias -> original contig name (identity if it was never aliased)."""
        return self.decode_map.get(str(alias), str(alias))

    def decode_variant_id(self, v_id: str) -> str:
        """Map a TRACE variant id ``[(LC)]{alias}_{pos}_{ref}_to_{alt}`` back.

        Parsed RIGHT-anchored, which is unambiguous because ``ref``/``alt`` are
        nucleotide strings and ``to`` is a literal. An id that does not have the
        expected shape is returned unchanged rather than mangled.
        """
        sp = str(v_id).split("_")
        if len(sp) < 5 or sp[-2] != "to":
            return str(v_id)
        head = "_".join(sp[:-4])
        prefix = ""
        if head.startswith("(LC)"):
            prefix, head = "(LC)", head[4:]
        return f"{prefix}{self.decode(head)}_" + "_".join(sp[-4:])

    @property
    def is_noop(self) -> bool:
        return not self.aliased


def alias_rows(
    rows: Iterable[MutableMapping[str, object]],
    alias: ContigAlias,
    key: str = "chrom",
) -> list[dict[str, object]]:
    """Return copies of `rows` with `key` replaced by its underscore-free alias.

    Copies rather than mutating, so a caller that still needs the original names (for
    a provenance record, or to decode the output) is not silently deprived of them.
    """
    out: list[dict[str, object]] = []
    for row in rows:
        new = dict(row)
        if key in new:
            new[key] = alias.encode(str(new[key]))
        out.append(new)
    return out
