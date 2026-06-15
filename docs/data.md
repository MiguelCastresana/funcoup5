# Data Notes

The full analyses require large external resources that are not committed to Git.

## Expected Local Layout

By default, `R/compare_network_properties.R` looks for:

```text
data/
├── benchmark/
│   ├── HumanNet-XN.tsv
│   └── gwas_catalog_v1.0-associations_e100_r2020-07-14.tsv
├── gold_standards/
│   ├── ppi_fc4
│   ├── ppi_FC5
│   ├── complex_fc4.txt
│   └── complex_FC5
├── networks/
│   ├── fc4.1_human
│   └── string_translated_0.8
└── genetranslations.csv
```

## Environment Variables

You can override those defaults:

```bash
export FUNCOUP_GOLD_DIR=/path/to/gold_standards
export FUNCOUP_BENCHMARK_DIR=/path/to/benchmark
export FUNCOUP_FC4_FILE=/path/to/fc4.1_human
export STRING_FILE=/path/to/string_translated_0.8
export GENE_TRANSLATIONS_FILE=/path/to/genetranslations.csv
```

## FunCoup Networks

Download current FunCoup networks from:

<https://funcoup.org/downloads/>

Keep downloaded network files in `data/networks/` or another ignored local directory.
