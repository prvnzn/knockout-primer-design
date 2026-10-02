# knockout-primers

Design the colony-PCR primers that **verify a Keio-collection gene deletion** in
*E. coli* — give it a list of gene names, get back a primer set that tells wild type
from knockout on a gel.

For every gene it designs two primer pairs, plus an optional third control:

| | where the primers sit | wild type | knockout |
|---|---|---|---|
| **Inner pair** | inside the ORF | band | **no product** |
| **Outer pair** | in the flanks either side | band | band, **shifted size** |
| **Cassette control** *(optional)* | inside the Kan&#8319; cassette + the matching flank | **no product** | band |

Together they give an unambiguous read: the inner pair proves the ORF is gone, the
outer pair proves you amplified the right locus (and shows the insertion as a size
shift), and the control proves that what replaced the ORF is the resistance cassette.

![How the primers bind in wild type and in the Keio deletion strain](docs/assets/assay.svg)

## The size shift, in one line

The Keio cassette (FRT–*aphA1*–FRT from pKD13) is a fixed **1323 bp**, and the flanks
are identical in both strains, so they cancel:

```
M_ko − M_wt  =  1323 − ORF length
```

The shift therefore depends **only on the ORF length** — not on where you put the outer
primers. A short ORF shifts strongly upward; an ORF close to 1323 bp barely moves. Genes
whose two bands would land within `outer.min_shift_bp` (150 bp by default) of each other
are flagged automatically, and for those you rely on the inner pair.

## Install

```bash
git clone https://github.com/prvnzn/knockout-primer-design.git
cd knockout-primer-design
pip install .
```

That gives you the `knockout-primers` command. If you would rather not install anything, just
`pip install -r requirements.txt` and run it as `python -m knockout_primers` instead — the two
are interchangeable everywhere below.

Needs Python 3.9+. The reference chromosome (~4.6 MB) is downloaded from NCBI the first
time you run it and cached in `cache/`; every run after that works offline.

## Use

Put your genes in a text file, one per line:

```
# proteases to verify
hslV
clpP
lon
degP
```

Then:

```bash
knockout-primers genes.txt --email you@example.com -o primers.csv
```

(or `python -m knockout_primers genes.txt --email you@example.com` without installing).
NCBI asks for an e-mail address with gene lookups — that is all `--email` is for.

```
  ok    hslV     -  ORF   531 bp  inner 522  outer 684->1476 (+792 bp), control 637 bp
  ok    clpP     +  ORF   624 bp  inner 612  outer 769->1468 (+699 bp), control 617 bp
  ok    lon      +  ORF  2355 bp  inner 2355 outer 2495->1463 (-1032 bp), control 557 bp
  check pepQ     +  ORF  1332 bp  inner 1329 outer 1478->1469 (-9 bp), control 568 bp
          - outer bands differ by only -9 bp (ORF 1332 bp vs cassette 1323 bp) — hard to
            tell apart on a gel, rely on the inner pair
```

`primers.csv` has four rows per gene (outer F/R, inner F/R) plus the control row, each with
the sequence, length, GC%, Tm, how many places in the genome it binds, both band sizes,
and any warning. It opens in Excel and pastes straight into an order form.

### Choosing your own criteria

Every threshold is yours to set. The common ones have flags:

```bash
knockout-primers genes.txt --tm-opt 62 --tm 58 68 --gc 45 65 --length 18 28
knockout-primers genes.txt --min-shift 200 --no-control
```

| flag | what it sets |
|---|---|
| `--tm-opt C` | melting temperature to aim for |
| `--tm MIN MAX` | Tm window to design within |
| `--tm-accept MIN MAX` | Tm range to accept before a primer is flagged |
| `--tm-pair-diff C` | largest Tm difference allowed within a pair |
| `--gc MIN MAX` | GC% window |
| `--length MIN MAX` | primer length, nt |
| `--cassette-bp BP` | cassette length, if you use something other than pKD13 |
| `--flank-bp BP` | how much genomic context to search for the outer pair |
| `--min-shift BP` | when to flag a gene as "bands too close to resolve" |
| `--no-control` | skip the cassette control |

Anything else in the config is reachable with `--set`:

```bash
knockout-primers genes.txt --set primers.max_poly_x=5 --set inner.min_coverage_frac=0.9
```

For a whole different set of criteria, take a copy and edit it — it is commented
throughout:

```bash
knockout-primers --write-config my_criteria.yaml
knockout-primers genes.txt --config my_criteria.yaml
```

### A CSV input, if you have more than names

Any of these columns, in any order (only `gene` is required):

```csv
gene,cds,outer_f,outer_r,inner_f,inner_r
hslV,,,,,
clpP,ATGTCATACAGCGGC...,,,,
degP,,GTTGACGCTGAAACTGG,CCAGTTTCAGCGTCAAC,ATGAAAAAAACCACATTAGCA,TGCTAATGTGGTTTTTTTCAT
```

- **`cds`** — locates the gene by aligning the sequence instead of by name. Use it when
  a name is ambiguous, or to pin a specific allele.
- **primer columns** — primers you already have. When all four are filled in, that gene is
  **scored and sized but never redesigned**, so a trusted, already-ordered set is safe.
  Leave them empty and the gene is designed fresh.

### From Python

```python
from knockout_primers.design import design_gene

r = design_gene("hslV")
r.outer_product_wt, r.outer_product_ko, r.outer_shift   # 684, 1476, +792
r.inner.fwd.seq, r.inner.fwd.tm                          # the oligo and its Tm
r.control.partner, r.control.product_bp                  # 'O_R', 637
r.warnings                                               # anything worth checking
```

## What it checks before handing you a primer

Each oligo is scored with [primer3](https://github.com/libnano/primer3-py) and flagged if
it misses any criterion: length, GC%, Tm and Tm match within the pair, hairpins,
self-dimers, cross-dimers, and runs of identical bases. On top of that:

- **Specificity** — every primer is matched against the whole chromosome; one that binds
  in more than one place is flagged, because it can give extra bands.
- **Placement** — the inner primers must genuinely lie inside the ORF, and the outer
  primers genuinely outside it. A primer in the wrong place silently breaks the assay, so
  it is checked rather than assumed.
- **Resolvable bands** — the outer shift must be large enough to see.

A gene marked `check` has a warning, not necessarily a problem. Read the `Notes` column
and decide — the tool does not quietly relax a threshold to make a gene pass.

## How the primers are chosen

Short version: the inner pair is pushed to the ORF's two ends so it spans as much of the
ORF as the quality criteria allow; the outer pair is anchored in the flanks and kept as
close to the ORF as possible, so both bands stay small and easy to run.
**[docs/METHODS.md](docs/METHODS.md)** sets out the whole procedure — why the assay is
built this way, exactly how each pair is placed and scored, how the band sizes and the
control are computed, and what the known limitations are.

## Beyond the Keio collection

Nothing is hard-coded to pKD13 or to *E. coli*. Any marker that replaces a complete ORF
works: set `cassette.length_bp` to its length and point `cassette.sequence_fasta` at its
sequence (needed only for the control), and set `control.sequence` to a primer inside your
cassette reading outward. Another organism needs its `reference.accession` and
`reference.organism` in the config.

## Tests

```bash
pip install pytest
pytest -m "not network"     # offline: config, parsing, assay geometry
pytest                      # also checks real genes against known-good values
```

The real-gene tests pin values produced by this assay's original implementation and
confirmed at the bench — `hslV`'s exact oligos, its 637 bp control band, and `pepQ`'s 9 bp
shift being flagged — so a regression in the geometry shows up as a failing test.

## Notes from the bench

- The Keio allele is modelled as a **clean, complete-ORF replacement** (Baba *et al.*
  2006). Real junctions carry a few scar bases, so sizes are accurate to within a handful
  of bp — fine for reading a gel, not a substitute for sequencing the junction.
- Gene coordinates come from **K-12 MG1655 (`NC_000913.3`)**. The Keio parent BW25113 is
  not in the NCBI Gene database but is sequence-identical at these loci.
- The control reports **both** the partner primer's label and its sequence. Forward/reverse
  conventions differ between labs, so set the tube up by sequence; if you get no band, try
  the other outer primer before suspecting the strain.
- Sequences are always handled and reported **5'→3'**.

## Citation

If this helped with a publication, please cite the tools it stands on — primer3
(Untergasser *et al.* 2012, *Nucleic Acids Res* 40:e115) and the Keio collection
(Baba *et al.* 2006, *Mol Syst Biol* 2:2006.0008).

## Licence

[MIT](LICENSE) — use it, change it, no warranty. Verify your own primers before ordering.
