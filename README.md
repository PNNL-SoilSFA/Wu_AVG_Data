# C12DO_MS
V-C12DO MS

This repository contains scripts and source data to reproduce analyses and figures for a global catechol 1,2-dioxygenase (C12DO) sequence embedding study, including (1) regression model scans of ESM2 embedding distance patterns, (2) a UMAP visualization highlighting a viral query (V-C12DO) and its nearest neighbors, and (3) residue conservation point plots and lineage level heatmaps anchored to aligned HMM columns.

## Contents

- `scripts/`
	- `c12do_analysis_pipeline.R`  
		Main analysis script, runs model scan, UMAP plot, and residue conservation plots.

- `data/`  
	All source input tables referenced by the script (listed below).

- `outputs/`  
	Generated tables and figures (created by the script).

## Source data

### Embedding distance and metadata tables

1. `c12do_esm2_t6_cpu_seed1.all_dist.tsv`  
	Required. Pairwise distance to the query embedding (or precomputed distances), with at least:
	- `seq_id`  
	- `cosine_distance`  
	- `length_count`

2. `c12do_esm2_t6_cpu_seed1.index_withTaxa_final.txt`  
	Required. Per sequence metadata, with at least:
	- `full_ID`  
	- `domain`  
	- `lineage_1`  
	- `lineage_2`  
	Optional:
	- `oligomer_simple` (used for labeling points on the UMAP)

### UMAP and neighbors tables

3. `c12do_esm2_t6_cpu_seed1.umap2d.tsv`  
	Required. 2D UMAP coordinates with:
	- `seq_id`  
	- `umap_x`  
	- `umap_y`

4. `c12do_esm2_t6_cpu_seed1.viral_neighbors.tsv`  
	Optional but recommended. Neighbor list for the query (V-C12DO).  
	The script auto-detects common column names for:
	- neighbor id column: `neighbor_id`, `neighbor`, `seq_id`, `full_ID`, `target_id`, `id`
	- distance column: `cosine_distance`, `distance`, `dist`, `cos_dist`

### Residue conservation and heatmap tables

These are keyed by `prefix`, currently:
- `prefix = c12do_matchcols_global`

5. `c12do_matchcols_global.matchcols.purity.tsv`  
	Required. Column occupancy or valid-count statistics per aligned HMM column, containing:
	- `match_col` (or one of: `Match_col`, `MatchCol`, `matchcol`, `HMM_col`)
	- occupancy column, one of: `occupancy`, `occupancy_valid`, `occ_valid`, `frac_valid`, `frac_non_gap`  
	or a valid-count column, one of: `n_valid`, `n_non_gap`, `n_obs`, `n`

6. `c12do_matchcols_global.matchcols.residue_composition.long.tsv`  
	Required. Long table of residue composition per `match_col`, containing:
	- `match_col` (or allowed aliases above)
	- amino acid column, one of: `aa`, `AA`, `residue`, `Residue`
	- value column, one of: `val`, `frac`, `fraction`, `freq`, `freq_valid`  
	Optional:
	- `n_valid` (per column), used to drop clearly unaligned columns.

7. `c12do_matchcols_global.keepcols.maxfreqvalid_ge_0.5.tsv`  
	Required. Subset of aligned columns used in the heatmap (must contain `match_col`).

8. `c12do_matchcols_global.perseq_residues.maxfreqvalid_ge_0.5.long.tsv`  
	Required. Per-sequence residue calls, long format, containing:
	- `seq_id`
	- `match_col`
	- `aa`

9. `c12do_keysites_rawtomatch_v6.keysites.per_seq.tsv`  
	Required. Sequence metadata table for the heatmap join, containing:
	- `seq_id`
	- `lineage_2`  
	Note: the script includes a small typo fix for a known lineage label variant.

10. `c12do_matchcols_global.keysites_from_alncols.tsv`  
	Optional but recommended. Key aligned HMM columns to highlight (must contain `match_col`).

## Outputs

The script writes the following files to the working directory unless you change paths.

### Model scan outputs
- `c12do_modelscan_length_domain_lineage.all_models.tsv`  
	All models evaluated, ranked by adjusted R2 then AIC.
- `c12do_modelscan_length_domain_lineage.best_model.tsv`  
	Best model row only.

### UMAP plot outputs
- `c12do_umap_domainShape_neighborFill_oligomerLabels.pdf`
- `c12do_umap_domainShape_neighborFill_oligomerLabels.png`

### Heatmap and intermediate tables
- `fig1b_vc_calls_table.tsv`
- `fig1b_heat_grid_table.tsv`
- `fig1b_heat_summarised_table.tsv`

## Requirements

- R (>= 4.1 recommended)
- R packages:
	- readr
	- dplyr
	- stats
	- ggplot2
	- ggrepel
	- scales
	- tidyr
	- stringr
	- patchwork
	- ggnewscale

Install packages in R:

install.packages(c(
	"readr","dplyr","ggplot2","ggrepel","scales",
	"tidyr","stringr","patchwork","ggnewscale"
))

## Running

Clone the repo:

git clone <this-repo>
cd <this-repo>


Put all input files under data/ (or update the file paths in the script).

Run the script:

Rscript scripts/c12do_analysis_pipeline.R

### Key parameters you may want to edit

In the script:

viral_id and vc_id
Query id used for the viral reference (default: V-C12DO).

k_neighbors
Number of neighbors for rays and highlighting.

min_n_per_level, max_levels_keep
Controls level lumping for regression predictors.

cutoff_frac
Residue conservation threshold used for alpha and labels.

lineage_n_min
Minimum number of called residues per lineage to include in heatmap.

Notes on identifiers

Sequence identifiers must match across files:

data/c12do_esm2_t6_cpu_seed1.all_dist.tsv uses seq_id

data/c12do_esm2_t6_cpu_seed1.index_withTaxa_final.txt uses full_ID
The script joins with by = c("seq_id" = "full_ID").

## Citation

If you use this repository, please cite the associated manuscript or dataset once available.

## License

Code: MIT License, data: CC BY 4.0.
