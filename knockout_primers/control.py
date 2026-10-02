"""The optional cassette (kanamycin) control.

The inner and outer pairs between them tell wild type from knockout, but neither
proves that what sits at the locus is *the cassette*. This control does: one
primer inside the resistance cassette, reading outward, paired with the outer
primer in the flank it points at.

    knockout:   flank ──► | cassette |          band
    wild type:  flank ──►   (no cassette)       no band

So it is the mirror image of the inner pair — a band only in the knockout. A
clean result is: inner band in wild type only, control band in the knockout only,
outer band in both at shifted sizes.

Which outer primer is the partner is pure geometry, worked out here from where
each primer actually binds the modelled knockout locus, rather than from a rule
of thumb. Both the primer and the expected size are reported, so there is no
ambiguity about which tube to set up.
"""
from __future__ import annotations

from dataclasses import dataclass

from .primers import Pair, Primer, score
from .reference import Locus, find_binding_site


@dataclass
class ControlResult:
    """The cassette-control PCR for one gene."""
    primer: Primer              # the cassette-internal primer (same for every gene)
    partner: str                # which outer primer pairs with it: 'O_F' or 'O_R'
    partner_seq: str
    product_bp: int             # expected band in the knockout
    note: str = ""

    @property
    def ok(self) -> bool:
        return self.product_bp > 0


def design_control(locus: Locus, outer: Pair | None, cassette: str,
                   cfg: dict) -> ControlResult | None:
    """Work out the cassette-control PCR for one gene.

    Returns ``None`` when the control is switched off, when no outer pair exists
    to partner with, or when the control primer does not bind the cassette.
    """
    ctrl_cfg = cfg.get("control") or {}
    if not ctrl_cfg.get("enabled", False) or outer is None:
        return None

    ctrl_seq = "".join((ctrl_cfg.get("sequence") or "").split()).upper()
    if not ctrl_seq:
        return None

    allele = locus.keio_allele(cassette)
    ctrl_site = find_binding_site(allele, ctrl_seq)
    if ctrl_site is None:
        return ControlResult(score(ctrl_seq, cfg), "", "", 0,
                             note="control primer does not match the cassette sequence")

    c_start, c_strand, c_len = ctrl_site
    # A primer extends away from its own 3' end: rightwards on the top strand,
    # leftwards on the bottom strand. Its partner must point back at it.
    for label, primer in (("O_F", outer.fwd), ("O_R", outer.rev)):
        site = find_binding_site(allele, primer.seq)
        if site is None:
            continue
        p_start, p_strand, p_len = site
        if p_strand == c_strand:
            continue                      # same strand: cannot be a convergent pair
        lo, hi = (p_start, c_start + c_len) if p_strand == "+" else (c_start, p_start + p_len)
        if hi <= lo:
            continue                      # divergent: they read away from each other
        return ControlResult(score(ctrl_seq, cfg), label, primer.seq, hi - lo)

    return ControlResult(score(ctrl_seq, cfg), "", "", 0,
                         note="neither outer primer forms a convergent pair with the "
                              "control primer — check the cassette orientation")
