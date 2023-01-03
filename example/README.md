## Examples

This folder contains example files for the remapping operations. They provide quick examples to showcase the functionality of this tool. Simply start up the Docker container and run these examples through using the `docker/run.sh` shell script (or `docker/map.R` directly).

For the examples below, note the following:
- ENSG00000121410 only maps to A1BG
- ENSG00000277796 maps to both CCL3L1 and CCL3L3
- ENSG00000261846 and ENSG00000197953 *both* map to AADACL2

#### `test_exp.tsv`:

This is a WebMeV-compatible expression matrix, which allows us to perform a MAD calculation for multi-mapping genes. Below, we map from Ensembl to gene symbol. The input is:

|gene|s1|s2|
|---|---|---|
|ENSG00000121410|1|2|
|ENSG00000277796|3|4|
|ENSG00000261846|5|6|
|ENSG00000197953|7|8|

and the output looks like:

|gene|s1|s2|
|---|---|---|
|A1BG|1|2|
|CCL3L1|3|4|
|CCL3L3|3|4|
|AADACL2|7|8|

Note that we get "duplicate" entries since ENSG00000277796 mapped to two gene symbols, CCL3L1 and CCL3L3.

Also (**important**) note that since both ENSG00000261846 and ENSG00000197953 map to AADACL2, we need to choose one or the other. To ensure that the mapping tool is compliant with other WebMeV analysis tools, our gene symbols *must* be unique. Hence, we can't have multipl AADACL2 entries. To remedy this, we apply the MAD calculation which selects the row corresponding to the ENSG00000197953 entry.

#### `test_ft.tsv`:

This is a WebMeV-compatible "feature table", which is used to hold information about (most often) genes. Generally, these are not strictly numeric tables, so we can't assume that we can perform a MAD calculation, or that it would be sensible (e.g. for a table of differential expression results with fold-change and p-values). Hence, our only recourse for creating WebMeV-compatible outputs is to take the first occurrence of any duplicates. The input looks like:

|gene|oncogene|pval|
|---|---|---|
|ENSG00000121410|Y|0.01|
|ENSG00000277796|N|0.05|
|ENSG00000261846|Y|0.001|
|ENSG00000197953|Y|0.85|

and the output looks like:

|gene|oncogene|pval|
|---|---|---|
|A1BG|Y|0.01|
|CCL3L1|N|0.05|
|CCL3L3|N|0.05|
|AADACL2|Y|0.85|

As before, we only have a single entry for gene symbol. *However*, note that due to a merge/join operation, the selected row is not strictly the first one that appears in the input file; the first row corresponding to AADACL2 had a p-value of 0.001, but the output table uses the row corresponding to the ENSG00000197953 Ensembl ID.