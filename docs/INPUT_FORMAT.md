# Input format

Every path in the config is resolved relative to `project_dir` unless it is
absolute. See `examples/` for working copies of every file described here.

## 1. Count matrix (required for qc_de, decomposition, entropy_pathways)

A CSV with one row per feature and one column per sample.

- **Column 1** is the feature identifier: a gene symbol (`A1BG`) or the
  DriverMap `LibID|SYMBOL` form (`2396|A1BG`), from which the symbol is taken.
- **Remaining columns** are samples, holding raw integer counts.
- Duplicate symbols are rejected unless `collapse_duplicate_symbols: true`, in
  which case their counts are summed.

```
feature,Parental_1,Parental_2,MetLung_1,MetLung_2
GENE0001,37,21,89,102
GENE0002,540,610,122,138
```

## 2. Sample sheet (required alongside the count matrix)

A CSV, one row per sample.

- `sample` (required): the analysis name for the sample.
- `group` (required): the biological state used to define contrasts.
- `context` (optional): free text such as `in_vitro` / `in_vivo`.
- `replicate` (optional): integer; filled in automatically if absent.
- `source_column` (optional): the exact column header in the count file that
  holds this sample. Use it when the count columns are not already named by
  `sample` (for example a raw instrument export with long descriptive
  headers). If omitted, the count columns must be named by `sample`.

```
sample,group,context,replicate
Parental_1,Parental,in_vitro,1
MetLung_1,MetLung,in_vivo,1
```

Contrasts in the config are written over the `group` values, as `"A - B"`,
read as "higher in A".

## 3. Dilution series (only for the crossreactivity module)

A CSV of raw counts, same feature-ID convention as the count matrix, with one
column per replicate of each mixture. Each mixture titrates target-species RNA
into a constant mass of host-species RNA. Column headers must contain the
`column_prefix` you list in the config; the module averages the replicates
sharing each prefix. You must include one pure-target mixture (`host_ng: 0`)
and one pure-host mixture (`target_ng: 0`).

```
feature,T50_H0_R1,T50_H0_R2,T10_H50_R1,T10_H50_R2,T0_H50_R1,T0_H50_R2
GENE0001,900,880,150,160,2,1
```

## 4. Optional resources (entropy_pathways)

- `ppi_edges`: a two-column, tab-separated edge list of gene symbols, for the
  network-entropy sub-module (SCENT SR and CCAT).
- `gene_sets`: a GMT file (the standard MSigDB format: set name, description,
  then tab-separated genes) for signature scoring.
- `collectri`: a CSV with columns `source,target,weight` for transcription-
  factor activity.

Any resource that is absent only disables its own sub-module; the rest still
run. PROGENy fetches its own model from the `progeny` R package and needs real
human gene symbols to match against.
