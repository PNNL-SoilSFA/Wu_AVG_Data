# ============================================================
# Cleaned script, redundant blocks removed, logic preserved
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stats)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(tidyr)
  library(stringr)
  library(patchwork)
  library(ggnewscale)
})

setwd("/Users/wuru978/Desktop/2026_projects/SFA_T1/second_round_dir/forVisual_dir/")

# ============================================================
# Part 1, model scan, effect sizes
# ============================================================

index_tsv <- "c12do_esm2_t6_cpu_seed1.index_withTaxa_final.txt"
dist_tsv <- "c12do_esm2_t6_cpu_seed1.all_dist.tsv"
out_prefix_model <- "c12do_modelscan_length_domain_lineage"

viral_id <- "V-C12DO"
min_n_per_level <- 10
max_levels_keep <- 30

lump_rare_levels <- function(x, min_n = 20, max_keep = 30, other_label = "Other") {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "Unknown"
  
  tab <- sort(table(x), decreasing = TRUE)
  keep_by_n <- names(tab)[tab >= min_n]
  keep_by_cap <- names(tab)[seq_len(min(length(tab), max_keep))]
  
  keep <- intersect(keep_by_n, keep_by_cap)
  if (length(keep) == 0) keep <- names(tab)[seq_len(min(length(tab), 2))]
  
  x[!(x %in% keep)] <- other_label
  factor(x)
}

fit_model <- function(df, formula_str) {
  m <- lm(as.formula(formula_str), data = df)
  s <- summary(m)
  list(r2 = s$r.squared, adj_r2 = s$adj.r.squared, aic = AIC(m))
}

df_idx <- read_tsv(index_tsv, show_col_types = FALSE)
df_dist <- read_tsv(dist_tsv, show_col_types = FALSE)

df <- df_dist %>%
  select(seq_id, cosine_distance, length_count) %>%
  left_join(
    df_idx %>% select(full_ID, domain, lineage_1, lineage_2),
    by = c("seq_id" = "full_ID")
  ) %>%
  filter(!is.na(cosine_distance), !is.na(length_count)) %>%
  mutate(is_query = (seq_id == viral_id))

df_stats <- df %>%
  filter(!is_query) %>%
  mutate(
    domain = lump_rare_levels(domain, min_n = min_n_per_level, max_keep = max_levels_keep),
    lineage_1 = lump_rare_levels(lineage_1, min_n = min_n_per_level, max_keep = max_levels_keep),
    lineage_2 = lump_rare_levels(lineage_2, min_n = min_n_per_level, max_keep = max_levels_keep)
  )

predictors <- c("length_count", "domain", "lineage_1", "lineage_2")

combo_list <- unlist(
  lapply(seq_along(predictors), function(i) combn(predictors, i, simplify = FALSE)),
  recursive = FALSE
)

base_r2 <- fit_model(df_stats, "cosine_distance ~ length_count")$r2

res_tbl <- bind_rows(lapply(combo_list, function(vars) {
  formula_str <- paste0("cosine_distance ~ ", paste(vars, collapse = " + "))
  f <- fit_model(df_stats, formula_str)
  tibble(
    model = formula_str,
    predictors = paste(vars, collapse = "+"),
    n_predictors = length(vars),
    r2 = f$r2,
    adj_r2 = f$adj_r2,
    aic = f$aic,
    delta_r2_vs_length = f$r2 - base_r2
  )
})) %>%
  arrange(desc(adj_r2), aic)

write_tsv(res_tbl, paste0(out_prefix_model, ".all_models.tsv"))

best <- res_tbl %>% slice(1)
best_formula <- best$model[[1]]
best_lm <- lm(as.formula(best_formula), data = df_stats)

write_tsv(best, paste0(out_prefix_model, ".best_model.tsv"))

a <- anova(best_lm)

term_rows <- setdiff(rownames(a), "Residuals")
ss_terms <- a[term_rows, "Sum Sq"]
names(ss_terms) <- term_rows
ss_res <- a["Residuals", "Sum Sq"]

eff <- tibble(
  term = names(ss_terms),
  ss = as.numeric(ss_terms),
  partial_eta2 = as.numeric(ss_terms / (ss_terms + ss_res)),
  pct_total_ss = as.numeric(ss_terms / sum(a[, "Sum Sq"]))
) %>%
  arrange(desc(partial_eta2))

print(best)
print(summary(best_lm))
print(eff)

# ============================================================
# Part 2, UMAP plot
# ============================================================

umap_tsv <- "c12do_esm2_t6_cpu_seed1.umap2d.tsv"
neighbors_tsv <- "c12do_esm2_t6_cpu_seed1.viral_neighbors.tsv"
out_prefix_umap <- "c12do_umap_domainShape_neighborFill_oligomerLabels"

k_neighbors <- 100

partial_eta2_lineage2 <- 0.49819409
partial_eta2_length <- 0.02375623
partial_eta2_domain <- 0.01222084

df_umap <- read_tsv(umap_tsv, show_col_types = FALSE)

df_plot <- df_umap %>%
  rename(umap1 = umap_x, umap2 = umap_y) %>%
  left_join(df_dist, by = "seq_id") %>%
  left_join(
    df_idx %>% select(full_ID, domain, lineage_2, oligomer_simple),
    by = c("seq_id" = "full_ID")
  ) %>%
  mutate(
    is_viral = (seq_id == viral_id),
    domain = factor(domain),
    lineage_2 = factor(lineage_2)
  )

stopifnot(sum(df_plot$is_viral) == 1)
viral_pt <- df_plot %>% filter(is_viral) %>% slice(1)

df_nei <- NULL
if (file.exists(neighbors_tsv)) {
  df_nei_raw <- read_tsv(neighbors_tsv, show_col_types = FALSE)
  
  id_candidates <- c("neighbor_id","neighbor","seq_id","full_ID","target_id","id")
  dist_candidates <- c("cosine_distance","distance","dist","cos_dist")
  
  id_col <- intersect(names(df_nei_raw), id_candidates)[1]
  dist_col <- intersect(names(df_nei_raw), dist_candidates)[1]
  
  if (is.na(id_col)) stop("Could not detect neighbor id column in neighbors_tsv.")
  if (is.na(dist_col)) stop("Could not detect distance column in neighbors_tsv.")
  
  df_nei <- df_nei_raw %>%
    rename(nei_id = all_of(id_col), nei_dist = all_of(dist_col)) %>%
    filter(!is.na(nei_id), !is.na(nei_dist)) %>%
    arrange(nei_dist) %>%
    slice_head(n = k_neighbors) %>%
    left_join(df_plot %>% select(seq_id, umap1, umap2, domain, lineage_2), by = c("nei_id" = "seq_id")) %>%
    filter(!is.na(umap1), !is.na(umap2))
} else {
  df_nei <- df_plot %>%
    filter(!is_viral, !is.na(cosine_distance)) %>%
    arrange(cosine_distance) %>%
    slice_head(n = k_neighbors) %>%
    transmute(nei_id = seq_id, nei_dist = cosine_distance, umap1, umap2, domain, lineage_2)
  message("neighbors_tsv not found, using smallest cosine_distance as neighbors.")
}

stopifnot(nrow(df_nei) > 0)

nei_ids <- unique(df_nei$nei_id)

df_plot <- df_plot %>%
  mutate(is_nei = seq_id %in% nei_ids)

df_seg <- df_nei %>%
  mutate(
    x = viral_pt$umap1,
    y = viral_pt$umap2,
    xend = umap1,
    yend = umap2,
    alpha_val = 1 - scales::rescale(nei_dist, to = c(0.05, 0.85))
  )

neighbor_lineages <- c(
  "Bacteria:Enterobacterales",
  "Bacteria:Kitasatosporales",
  "Bacteria:Mycobacteriales",
  "Bacteria:Pseudomonadales",
  "Bacteria:Pseudonocardiales",
  "Bacteria:Rhodospirillales",
  "Bacteria:Sphingomonadales",
  "Bacteria:Lysobacterales",
  "Bacteria:Neisseriales",
  "Bacteria:Micrococcales"
)

cols_nei <- c(
  "Bacteria:Enterobacterales" = "#D55E00",
  "Bacteria:Kitasatosporales" = "#56B4E9",
  "Bacteria:Mycobacteriales" = "#999999",
  "Bacteria:Pseudomonadales" = "#E69F00",
  "Bacteria:Pseudonocardiales" = "#009E73",
  "Bacteria:Rhodospirillales" = "#0072B2",
  "Bacteria:Sphingomonadales" = "#CC79A7",
  "Bacteria:Lysobacterales" = "#756bb1",
  "Bacteria:Neisseriales" = "#ffeda0",
  "Bacteria:Micrococcales" = "#A6761D"
)

df_plot <- df_plot %>%
  mutate(
    lineage_2_nei = ifelse(is_nei, as.character(lineage_2), "non_neighbor"),
    lineage_2_nei = factor(lineage_2_nei, levels = c(neighbor_lineages, "non_neighbor"))
  )

cols_fill <- c(cols_nei, non_neighbor = "white")

df_oligo <- df_plot %>%
  filter(!is.na(oligomer_simple), oligomer_simple != "unknown") %>%
  distinct(seq_id, .keep_all = TRUE) %>%
  mutate(fill_group = ifelse(is_nei, as.character(lineage_2), "non_neighbor"))

ann_text <- paste0(
  "Variance explained by\n",
  " lineage of the 10,950 C12DOs = ", sprintf("%.3f", partial_eta2_lineage2), "\n",
  " sequence length = ", sprintf("%.3f", partial_eta2_length), "\n",
  " domain category = ", sprintf("%.3f", partial_eta2_domain)
)

shape_map <- c(
  "Bacteria" = 21,
  "Eukaryota:Arthropoda" = 22,
  "Eukaryota:Fungi" = 24,
  "Viruses" = 23
)

p_umap <- ggplot() +
  geom_point(
    data = df_plot %>% filter(!is_viral, !is_nei),
    aes(x = umap1, y = umap2, shape = domain),
    color = "grey50",
    fill = NA,
    stroke = 0.35,
    size = 0.7,
    alpha = 0.30
  ) +
  geom_segment(
    data = df_seg,
    aes(x = x, y = y, xend = xend, yend = yend, alpha = alpha_val),
    color = "grey35",
    linewidth = 0.25
  ) +
  scale_alpha_continuous(guide = "none") +
  geom_point(
    data = df_plot %>% filter(domain %in% c("Eukaryota:Arthropoda","Eukaryota:Fungi"), !is_viral, !is_nei),
    aes(x = umap1, y = umap2, shape = domain),
    color = "grey35",
    fill = NA,
    stroke = 0.5,
    size = 1,
    alpha = 1
  ) +
  geom_point(
    data = df_plot %>% filter(!is_viral, is_nei),
    aes(x = umap1, y = umap2, shape = domain, fill = lineage_2_nei),
    color = "black",
    stroke = 0.25,
    size = 1.4,
    alpha = 0.95
  ) +
  scale_fill_manual(values = cols_fill, drop = FALSE) +
  geom_point(
    data = df_plot %>% filter(is_viral),
    aes(x = umap1, y = umap2, shape = domain),
    color = "black",
    fill = "red",
    stroke = 0.6,
    size = 3,
    alpha = 1
  ) +
  ggrepel::geom_text_repel(
    data = df_plot %>% filter(is_viral),
    aes(x = umap1 + 0.8, y = umap2 + 0.4, label = seq_id),
    size = 3,
    box.padding = 0.4,
    point.padding = 0.1
  ) +
  geom_point(
    data = df_oligo,
    aes(x = umap1, y = umap2, fill = fill_group),
    shape = 21,
    size = 1.2,
    stroke = 0.9,
    color = "black",
    alpha = 1
  ) +
  ggrepel::geom_text_repel(
    data = df_oligo,
    aes(x = umap1, y = umap2, label = oligomer_simple),
    size = 3,
    seed = 1,
    nudge_y = -0.35,
    box.padding = 0.4,
    point.padding = 0.1,
    force_pull = 0.2,
    min.segment.length = 0,
    segment.size = 0.3,
    max.overlaps = Inf
  ) +
  scale_shape_manual(values = shape_map) +
  theme_classic() +
  labs(x = "UMAP1", y = "UMAP2", shape = "Domain", fill = "Neighbor lineage") +
  guides(
    fill = guide_legend(
      title = "Neighbor lineage",
      override.aes = list(shape = 21, size = 2.8, alpha = 1),
      ncol = 1
    ),
    shape = guide_legend(
      title = "Domain",
      override.aes = list(fill = "white", color = "black", size = 2.8, alpha = 1)
    )
  ) +
  annotate(
    "label",
    x = -9.5, y = 18,
    label = ann_text,
    hjust = 0,
    vjust = 1,
    size = 2.5,
    linewidth = 0.2
  ) +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

print(p_umap)
ggsave(paste0(out_prefix_umap, ".pdf"), p_umap, width = 10.5, height = 6.5)
ggsave(paste0(out_prefix_umap, ".png"), p_umap, width = 10.5, height = 6.5, dpi = 300)

# ============================================================
# Part 3, residue conservation plots, points plus heatmap
# ============================================================

prefix <- "c12do_matchcols_global"

purity_tsv <- paste0(prefix, ".matchcols.purity.tsv")
comp_tsv <- paste0(prefix, ".matchcols.residue_composition.long.tsv")
keymap_tsv <- paste0(prefix, ".keysites_from_alncols.tsv")

meta_all_tsv <- "c12do_keysites_rawtomatch_v6.keysites.per_seq.tsv"
perseq_tsv <- paste0(prefix, ".perseq_residues.maxfreqvalid_ge_0.5.long.tsv")
keepcols_tsv <- paste0(prefix, ".keepcols.maxfreqvalid_ge_0.5.tsv")

vc_id <- "V-C12DO"

occ_cutoff <- 0
full_frac_threshold <- 0.999999
cutoff_frac <- 0.50
alpha_strong <- 1.00
alpha_faint <- 0.07

segment_width <- 110
pretty_n_breaks <- 15
point_size <- 0.60
point_jitter_width <- 0.08

label_keysites_always <- TRUE
label_y <- 1.04
label_text_size <- 2.4

lineage_n_min <- 1
heat_x_label_every <- 2

aa_levels_fixed <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
aa20_colors <- c(
  "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
  "#8B8000","#A65628","#F781BF","#999999","#088F8F",
  "#FC8D62","#0000FF","#E78AC3","#A6D854","#FFD92F",
  "#E5C494","#B3B3B3","#1B9E77","#D95F02","#7570B3"
)
pal <- setNames(aa20_colors, aa_levels_fixed)

fill_low <- "#89CFF0"
fill_mid <- "white"
fill_high <- "#FA8072"

col_vc_missing <- "grey92"
col_lin_missing <- "grey35"

safe_read_tsv <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  read_tsv(path, show_col_types = FALSE, progress = FALSE)
}

standardize_match_col <- function(df, colname_guess = c("match_col","Match_col","MatchCol","matchcol","HMM_col")) {
  hit <- intersect(colname_guess, names(df))
  if (length(hit) == 0) stop("Could not find match_col in: ", paste(names(df), collapse = ", "))
  df %>% mutate(match_col = as.integer(.data[[hit[1]]]))
}

fix_lineage_typos <- function(x) {
  x <- str_trim(as.character(x))
  x[x == "Bacteria:unclassfiedpseudomonadati"] <- "Bacteria:unclassifiedpseudomonadati"
  x
}

df_purity <- safe_read_tsv(purity_tsv) %>% standardize_match_col()
df_comp_long_raw <- safe_read_tsv(comp_tsv) %>% standardize_match_col()

df_keymap <- if (file.exists(keymap_tsv)) safe_read_tsv(keymap_tsv) %>% standardize_match_col() else NULL

aa_col <- intersect(names(df_comp_long_raw), c("aa","AA","residue","Residue"))[1]
val_col <- intersect(names(df_comp_long_raw), c("val","frac","fraction","freq","freq_valid"))[1]
if (is.na(aa_col) || is.na(val_col)) stop("Could not find aa and value columns in comp table.")

df_comp_long <- df_comp_long_raw %>%
  mutate(
    aa = toupper(as.character(.data[[aa_col]])),
    val = as.numeric(.data[[val_col]])
  ) %>%
  filter(!is.na(match_col), !is.na(aa), is.finite(val)) %>%
  filter(!(aa %in% c("OTHER","Other"))) %>%
  filter(aa %in% aa_levels_fixed)

df_sum <- df_comp_long %>% group_by(match_col) %>% summarise(s = sum(val, na.rm = TRUE), .groups = "drop")
need_norm <- !(mean(abs(df_sum$s - 1) < 1e-3, na.rm = TRUE) > 0.95)
if (need_norm) {
  df_comp_long <- df_comp_long %>%
    group_by(match_col) %>%
    mutate(val = val / sum(val, na.rm = TRUE)) %>%
    ungroup()
}

occ_col <- intersect(names(df_purity), c("occupancy","occupancy_valid","occ_valid","frac_valid","frac_non_gap"))[1]
n_col <- intersect(names(df_purity), c("n_valid","n_non_gap","n_obs","n"))[1]

df_allowed <- if (!is.na(occ_col)) {
  df_purity %>%
    mutate(occ = as.numeric(.data[[occ_col]])) %>%
    filter(is.finite(occ), occ >= occ_cutoff) %>%
    distinct(match_col)
} else if (!is.na(n_col)) {
  df_purity %>%
    mutate(nv = as.numeric(.data[[n_col]])) %>%
    filter(is.finite(nv), nv > 0) %>%
    distinct(match_col)
} else {
  stop("Could not find occupancy or n_valid in purity table.")
}

df_comp_long <- df_comp_long %>% semi_join(df_allowed, by = "match_col")

df_pos <- df_comp_long %>%
  distinct(match_col) %>%
  arrange(match_col) %>%
  mutate(
    plot_pos = row_number(),
    segment = (plot_pos - 1) %/% segment_width + 1
  )

df_pos2 <- df_pos %>% distinct(segment, plot_pos, match_col)

df_ticks <- df_pos2 %>%
  group_by(segment) %>%
  arrange(plot_pos) %>%
  mutate(step = max(1L, ceiling(n() / pretty_n_breaks))) %>%
  filter((row_number() - 1) %% first(step) == 0) %>%
  ungroup() %>%
  select(plot_pos, match_col)

df_plot_full <- df_pos2 %>%
  crossing(aa = aa_levels_fixed) %>%
  left_join(df_comp_long %>% select(match_col, aa, val), by = c("match_col","aa")) %>%
  mutate(
    val = replace_na(val, 0),
    aa = factor(aa, levels = aa_levels_fixed),
    alpha_val = ifelse(val >= cutoff_frac, alpha_strong, alpha_faint)
  )

df_full <- df_plot_full %>%
  group_by(match_col) %>%
  summarise(max_frac = max(val, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(max_frac), max_frac >= full_frac_threshold) %>%
  select(match_col)

df_key <- NULL
if (!is.null(df_keymap) && nrow(df_keymap) > 0) {
  df_key <- df_keymap %>% distinct(match_col) %>% semi_join(df_allowed, by = "match_col")
  if (nrow(df_key) > 4) df_key <- df_key %>% arrange(match_col) %>% slice(1:4)
}

df_labels <- df_plot_full %>%
  group_by(segment, plot_pos, match_col) %>%
  slice_max(val, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  filter(is.finite(val), val >= cutoff_frac) %>%
  mutate(aa = as.character(aa))

if (isTRUE(label_keysites_always) && !is.null(df_key) && nrow(df_key) > 0) {
  df_key2 <- df_key %>% inner_join(df_pos2, by = "match_col")
  df_key_labels <- df_plot_full %>%
    inner_join(df_key2 %>% select(segment, plot_pos, match_col), by = c("segment","plot_pos","match_col")) %>%
    group_by(segment, plot_pos, match_col) %>%
    slice_max(val, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(aa = as.character(aa))
  df_labels <- bind_rows(df_labels, df_key_labels) %>% distinct(segment, plot_pos, match_col, aa, .keep_all = TRUE)
}

p_comp <- ggplot(df_plot_full, aes(x = plot_pos, y = val, color = aa)) +
  geom_point(aes(alpha = alpha_val), size = point_size, position = position_jitter(width = point_jitter_width, height = 0)) +
  facet_wrap(~ segment, scales = "free_x", ncol = 1) +
  scale_y_continuous(limits = c(0, 1.05), labels = percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = df_ticks$plot_pos, labels = df_ticks$match_col) +
  geom_text(
    data = df_labels,
    aes(x = plot_pos, y = label_y, label = aa, color = aa),
    size = label_text_size,
    show.legend = FALSE
  ) +
  scale_color_manual(values = pal) +
  scale_alpha_identity(guide = "none") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = "Aligned HMM column position", y = "Residue fraction", color = NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.35, "cm"),
    strip.text = element_blank(),
    strip.background = element_blank()
  )

if (nrow(df_full) > 0) {
  p_comp <- p_comp + geom_vline(
    data = df_full %>% inner_join(df_pos2, by = "match_col"),
    aes(xintercept = plot_pos),
    linetype = "dashed",
    linewidth = 0.25,
    color = "grey30"
  )
}

if (!is.null(df_key) && nrow(df_key) > 0) {
  p_comp <- p_comp + geom_vline(
    data = df_key %>% inner_join(df_pos2, by = "match_col"),
    aes(xintercept = plot_pos),
    linetype = "dashed",
    linewidth = 0.4,
    color = "red3"
  )
}

df_keepcols <- safe_read_tsv(keepcols_tsv) %>%
  mutate(match_col = as.integer(match_col)) %>%
  arrange(match_col)

x_levels <- df_keepcols$match_col

df_perseq <- safe_read_tsv(perseq_tsv) %>%
  transmute(
    seq_id = as.character(seq_id),
    match_col = as.integer(match_col),
    aa = toupper(as.character(aa))
  ) %>%
  filter(match_col %in% x_levels)

df_meta <- safe_read_tsv(meta_all_tsv) %>%
  transmute(seq_id = as.character(seq_id), lineage_2 = fix_lineage_typos(lineage_2)) %>%
  filter(!is.na(lineage_2), lineage_2 != "")

df_vc_calls <- df_perseq %>% filter(seq_id == vc_id) %>% distinct(match_col, vc_aa = aa)

df_vc_all <- tibble(match_col = x_levels) %>%
  left_join(df_vc_calls, by = "match_col") %>%
  mutate(
    vc_missing = is.na(vc_aa),
    vc_aa_plot = ifelse(vc_missing, "-", vc_aa),
    x = factor(match_col, levels = x_levels)
  )

if (all(df_vc_all$vc_missing)) stop("V-C12DO has no residue calls, check vc_id or perseq_tsv.")

df_perseq2 <- df_perseq %>%
  inner_join(df_meta, by = "seq_id") %>%
  inner_join(df_vc_all %>% filter(!vc_missing) %>% select(match_col, vc_aa), by = "match_col")

df_heat <- df_perseq2 %>%
  group_by(lineage_2, match_col, vc_aa) %>%
  summarise(n_called = n(), match_frac = mean(aa == vc_aa), .groups = "drop")

df_lineage_n <- df_heat %>%
  group_by(lineage_2) %>%
  summarise(n_total_called = sum(n_called), .groups = "drop")

keep_lineages <- df_lineage_n %>%
  filter(n_total_called >= lineage_n_min) %>%
  pull(lineage_2)

df_grid <- expand.grid(
  lineage_2 = keep_lineages,
  match_col = df_vc_all %>% filter(!vc_missing) %>% pull(match_col),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  left_join(df_heat %>% select(lineage_2, match_col, match_frac, n_called), by = c("lineage_2","match_col")) %>%
  mutate(
    n_called = ifelse(is.na(n_called), 0L, n_called),
    match_frac = ifelse(n_called == 0L, NA_real_, match_frac),
    x = factor(match_col, levels = x_levels)
  )

y_levels <- df_grid %>%
  group_by(lineage_2) %>%
  summarise(mean_match = mean(match_frac, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_match)) %>%
  pull(lineage_2)

df_heat_plot <- df_grid %>%
  filter(!is.na(match_frac)) %>%
  mutate(lineage_2 = factor(lineage_2, levels = y_levels))

df_missing_lin <- df_grid %>%
  filter(is.na(match_frac)) %>%
  transmute(
    lineage_2 = factor(lineage_2, levels = y_levels),
    x = x,
    miss_type = factor("Lineage gap", levels = c("V-C12DO gap","Lineage gap"))
  )

df_missing_vc_tile <- df_vc_all %>%
  filter(vc_missing) %>%
  select(x) %>%
  distinct() %>%
  crossing(lineage_2 = factor(y_levels, levels = y_levels)) %>%
  mutate(miss_type = factor("V-C12DO gap", levels = c("V-C12DO gap","Lineage gap")))

df_top <- df_vc_all %>%
  mutate(y = 1, lab_col = ifelse(vc_missing, "grey60", "black"))

x_right <- nlevels(df_top$x) + 0.6

p_top <- ggplot(df_top, aes(x = x, y = y)) +
  geom_text(aes(label = vc_aa_plot), size = 2, color = df_top$lab_col) +
  annotate("text", x = x_right, y = 1, label = vc_id, hjust = 0, size = 2.8) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(t = 4, r = 50, b = -10, l = 10))

show_idx <- seq_along(x_levels)[(seq_along(x_levels) - 1) %% heat_x_label_every == 0]
x_breaks <- x_levels[show_idx]

p_heat <- ggplot() +
  geom_tile(data = df_missing_vc_tile, aes(x = x, y = lineage_2, fill = miss_type)) +
  geom_tile(data = df_missing_lin, aes(x = x, y = lineage_2, fill = miss_type)) +
  scale_fill_manual(
    values = c("V-C12DO gap" = col_vc_missing, "Lineage gap" = col_lin_missing),
    name = "Gap position",
    drop = FALSE
  ) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = df_heat_plot, aes(x = x, y = lineage_2, fill = match_frac)) +
  scale_fill_gradient2(
    low = fill_low,
    mid = fill_mid,
    high = fill_high,
    midpoint = 0.5,
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = percent_format(accuracy = 1),
    name = "Percent of lineage sequences\nmatched to V-C12DO"
  ) +
  scale_x_discrete(drop = FALSE, breaks = as.character(x_breaks), labels = as.character(x_breaks)) +
  labs(x = "Aligned HMM column position", y = "Lineage") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 6.5, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(t = -20, r = 10, b = 6, l = 10)
  )

fig_heat <- p_top / p_heat + plot_layout(heights = c(0.06, 1))
fig_combo1 <- p_umap / fig_heat + plot_layout(heights = c(0.6, 1))

print(fig_combo1)

write_tsv(df_vc_all, "fig1b_vc_calls_table.tsv")
write_tsv(df_grid, "fig1b_heat_grid_table.tsv")
write_tsv(df_heat, "fig1b_heat_summarised_table.tsv")
