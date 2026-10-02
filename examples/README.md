# Example inputs

**`genes.txt`** — the usual case: one gene name per line.

```bash
knockout-primers examples/genes.txt --email you@example.com -o primers.csv
```

**`genes_with_options.csv`** — the CSV form, when you have more than names. Only `gene` is
required; the other columns may be left empty.

- `cds` — locate the gene by aligning this sequence instead of by name.
- `outer_f` / `outer_r` / `inner_f` / `inner_r` — primers you already have. Fill in all
  four and that gene is scored and sized but **never redesigned**.

```bash
knockout-primers examples/genes_with_options.csv --email you@example.com -o primers.csv
```
