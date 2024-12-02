#############################
### Install/load packages ###
#############################

#install.packages("librarian")
librarian::shelf(tidyverse, factoextra, data.table, mpra, reshape2,
                 ggtext, UpSetR, patchwork, GGally,
                 EnsDb.Hsapiens.v79, locuszoomr, ensembldb, Biostrings,
                 ggrepel, gkmSVM, rstatix, ggpubr, scales, BSgenome.Hsapiens.UCSC.hg38,
                 chromVARmotifs, ggh4x, memes)

#####################################
### Functions for processing data ###
#####################################

filtActive <- function(res){
  res %>% dplyr::filter(adj.P.Val < 0.05) %>% dplyr::select(new_refname)
}

filtNonActive <- function(res){
  res %>% dplyr::filter(adj.P.Val > 0.05) %>% dplyr::select(new_refname)
}

processPairs = function(df) {
  df %>%
    extract(new_refname, into = c("refname", "prom"), "(.*)_([^_]+)$") %>%
    separate(refname, into = c("rsid", "chr", "pos", "ref", "alt",
                               "allele", "class", "site"),
             sep = "_", remove = FALSE) %>%
    mutate(class = case_when(ref == "NA" & alt == "NA" & allele == "NA" ~ "neg_ctrl",
                             alt == "plus" | alt == "minus" ~ "cage_strand",
                             class == "left" & is.na(site) ~ "standard_left",
                             class == "right" & is.na(site) ~ "standard_right",
                             site == "left" & class == "del" ~ "deletion_left",
                             site == "right" & class == "del" ~ "deletion_right",
                             site == "left" & class == "shf" ~ "shuffle_left",
                             site == "right" & class == "shf" ~ "shuffle_right",
                             class == "del" & is.na(site) ~ "deletion",
                             class == "shf" & is.na(site) ~ "shuffle",
                             TRUE ~ "standard"),
           allele = case_when(alt == "minus" ~ "A",
                              alt == "plus" ~ "R",
                              TRUE ~ allele),
           rsid = case_when(class == "cage_strand" ~ pos,
                            TRUE ~ rsid),
           rsid = paste(rsid, class, sep = "_"),
           new_refname = paste(rsid, prom, sep = "_")) %>%
    dplyr::select(new_refname, refname, prom, rsid, allele, chr, pos, class, site)
}

##################
### Set up env ###
##################

# set up output directories
work_dir <- "/path/to/workdir/mpra"

# in data dir - should have mpra counts and dictionary
in_dir <- "/path/to/data"

# set fig dir
fig_dir <- file.path(work_dir, "figs")
dir.create(fig_dir) # only run once

# out data dir
out_dir <- file.path(work_dir, "output")
dir.create(out_dir) # only run once

############################################
### Load dictionaries and raw count data ###
############################################

# Load in pairing dictionary
full_dict = fread(file.path(in_dir, "tovar-nishino_library_barcode_pairing.txt"))

# Load in raw counts
raw_cts = fread(file.path(in_dir, "joint_matrix.txt"))
colnames(raw_cts)[c(2:34)] <- c(paste0(c("dna", "rna"), rep(1:16, each = 2)),"nrt")

####################
### Process data ###
####################

# Merge raw counts with full dictionary
merged_raw_cts <- merge(full_dict, raw_cts, by = "barcode")

# Add counts -- aka number of barcodes per oligo
merged_raw_cts <- merged_raw_cts %>%
  group_by(refname, prom) %>%
  add_count(refname)

# Filter raw counts for any oligos with 10+ unique barcodes
filt_cts <- merged_raw_cts %>%
  dplyr::filter(n > 10) %>%
  group_by(refname, prom) %>%
  reframe(refname = refname,
          prom = prom,
          across(c(dna1:rna16), sum)) %>%
  distinct(refname, prom, .keep_all = TRUE)

write.table(filt_cts, file.path(out_dir, "filt_cts.tsv"),
            sep = "\t", quote = F, row.names = F)

# Calculate TPM values for all oligos (library size-normalized)
filt_tpm <- filt_cts %>% 
  pivot_longer(starts_with(c("rna","dna")), 
               names_to = "sample", 
               values_to = "counts") %>%
  dplyr::select(refname, prom, sample, counts) %>%
  group_by(sample) %>%
  mutate(tpm = (counts / (sum(counts)) * 1e6)) %>%
  dplyr::select(refname, prom, sample, tpm) %>%
  pivot_wider(names_from = sample, values_from = tpm)

# Calculate log2(RNA/DNA) for later plotting purposes
plot_tpm <- filt_tpm %>%
  mutate(ratio1 = log2(rna1/dna1),
         ratio2 = log2(rna2/dna2),
         ratio3 = log2(rna3/dna3),
         ratio4 = log2(rna4/dna4),
         ratio5 = log2(rna5/dna5),
         ratio6 = log2(rna6/dna6),
         ratio7 = log2(rna7/dna7),
         ratio8 = log2(rna8/dna8),
         ratio9 = log2(rna9/dna9),
         ratio10 = log2(rna10/dna10),
         ratio11 = log2(rna11/dna11),
         ratio12 = log2(rna12/dna12),
         ratio13 = log2(rna13/dna13),
         ratio14 = log2(rna14/dna14),
         ratio15 = log2(rna15/dna15),
         ratio16 = log2(rna16/dna16))

#################
#### QC plots ###
#################

# Check pairwise correlation
pairplot <- ggpairs(log10(plot_tpm[,c(35:50)]+1),
                    lower = list(continuous = wrap("points", alpha = 0.1, size=0.5), 
                                 combo = wrap("dot", alpha = 0.01, size=0.6)))

# PCA
pca_df <- plot_tpm[complete.cases(plot_tpm[,c(35:50)]),]
pca_df <- pca_df[is.finite(rowSums(pca_df[,c(3:50)])),]

pca_res <- prcomp(t(pca_df[,c(35:50)]))  # Transpose to have samples in rows
# Create a data frame with PCA results
pca_res_df <- data.frame(
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2],
  condition = factor(c(rep("undiff",4),rep("diff",4),rep("aicar",4),rep("palm",4)))
)


# Calculate the proportion of variance explained
pct_var_pca <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), digits = 2)

# Plot
pca_res_df$label <- c(NA, "Undifferentiated\n+ Basal", rep(NA, 4), "Differentiated\n+ Basal",
                      rep(NA, 3), "Differentiated\n+ AICAR", rep(NA, 3), "Differentiated\n+ Palmitate", NA)
pca_res_df$color <- c(rep("#ff595e",4),rep("#ffca3a",4),rep("#8ac926",4), rep("#1982c4", 4))
pca_res_df$status <- c(rep("undiff",4),rep("diff",12))
mpra_pca <- ggplot(pca_res_df, aes(PC1, PC2, color = condition, label=label, shape=status)) +
  geom_point(size=3, fill = pca_res_df$color, color = "black") +
  scale_shape_manual(values = c(21, 24)) +
  xlab(paste0("PC1: ",pct_var_pca[1],"% variance")) +
  ylab(paste0("PC2: ",pct_var_pca[2],"% variance")) + 
  scale_color_manual(values = c("#ff595e","#ffca3a","#8ac926","#1982c4")) +
  geom_text_repel(min.segment.length = 0, nudge_y = 0.5, family = "Helvetica",
                  box.padding = 1.5, size = 4, lineheight = 1,
                  data = pca_res_df[c(1:4),], color = "black") + 
  geom_text_repel(min.segment.length = 0, family = "Helvetica",
                  box.padding = 1.5, nudge_y = -0.1, nudge_x = -0.1, size = 4, lineheight = 1,
                  data = pca_res_df[c(5:8),], color = "black") +
  geom_text_repel(min.segment.length = 0, family = "Helvetica",
                  box.padding = 1.5, nudge_y = -0.1, nudge_x = -0.1, size = 4, lineheight = 1,
                  data = pca_res_df[c(9:12),], color = "black") +
  geom_text_repel(min.segment.length = 0, family = "Helvetica",
                  box.padding = 1.5, nudge_x = -0.05, nudge_y = 0.05, size = 4, lineheight = 1,
                  data = pca_res_df[c(13:16),], color = "black") +
  theme_bw(base_family = "Helvetica", base_size = 14) + theme(legend.position="none")

ggsave(mpra_pca, filename = file.path(fig_dir, "mpra-pca.png"),
       units = "in", dpi = 600, width = 4, height = 3.2, device = ragg::agg_png())

# do pca for promoters separately
## myb
plot_tpm_myb <- plot_tpm %>% dplyr::filter(prom == "MYBPC2")
pca_df_myb <- plot_tpm_myb[complete.cases(plot_tpm_myb[,c(35:50)]),]
pca_df_myb <- pca_df_myb[is.finite(rowSums(pca_df_myb[,c(3:50)])),]

pca_res_myb <- prcomp(t(pca_df_myb[,c(35:50)]))  # Transpose to have samples in rows
# Create a data frame with PCA results
pca_res_myb_df <- data.frame(
  PC1 = pca_res_myb$x[,1],
  PC2 = pca_res_myb$x[,2],
  condition = factor(c(rep("undiff",4),rep("diff",4),rep("aicar",4),rep("palm",4)))
)


# Calculate the proportion of variance explained
pct_var_myb_pca <- round(100 * (pca_res_myb$sdev^2 / sum(pca_res_myb$sdev^2)), digits = 2)

# Plot
pca_res_myb_df$label <- c(NA, "Undifferentiated\n+ Basal", rep(NA, 4), "Differentiated\n+ Basal",
                          rep(NA, 3), "Differentiated\n+ AICAR", rep(NA, 3), "Differentiated\n+ Palmitate", NA)
pca_res_myb_df$color <- c(rep("#ff595e",4),rep("#ffca3a",4),rep("#8ac926",4), rep("#1982c4", 4))
pca_res_myb_df$status <- c(rep("undiff",4),rep("diff",12))
mpra_myb_pca <- ggplot(pca_res_myb_df, aes(PC1, PC2, color = condition, label=label, shape=status)) +
  geom_point(size=3, fill = pca_res_myb_df$color, color = "black") +
  scale_shape_manual(values = c(21, 24)) +
  xlab(paste0("PC1: ",pct_var_myb_pca[1],"% variance")) +
  ylab(paste0("PC2: ",pct_var_myb_pca[2],"% variance")) + 
  scale_color_manual(values = c("#ff595e","#ffca3a","#8ac926","#1982c4")) +
  geom_text_repel(min.segment.length = 0, family = "Helvetica",
                  box.padding = 1.5, size = 5, lineheight = 1,
                  data = pca_res_myb_df[c(1:4),], color = "black") + 
  geom_text_repel(min.segment.length = 0, family = "Helvetica",
                  box.padding = 1.5, nudge_y = -0.05, size = 5, lineheight = 1,
                  data = pca_res_myb_df[c(5:8),], color = "black") +
  geom_text_repel(min.segment.length = 0, family = "Helvetica",
                  box.padding = 1.5, nudge_y = -0.4, size = 5, lineheight = 1,
                  data = pca_res_myb_df[c(9:12),], color = "black") +
  geom_text_repel(min.segment.length = 0, family = "Helvetica",
                  box.padding = 1.5, nudge_x = -0.05, size = 5, lineheight = 1,
                  data = pca_res_myb_df[c(13:16),], color = "black") +
  theme_bw(base_family = "Helvetica", base_size = 16) + theme(aspect.ratio = 1, legend.position="none")

ggsave(mpra_myb_pca, filename = file.path(fig_dir, "mpra-myb-pca.png"),
       units = "in", dpi = 600, width = 5, height = 5, device = ragg::agg_png())

# scp
plot_tpm_scp <- plot_tpm %>%
  dplyr::filter(prom == "SCP1")
pca_df_scp <- plot_tpm_scp[complete.cases(plot_tpm_scp[,c(35:50)]),]
pca_df_scp <- pca_df_scp[is.finite(rowSums(pca_df_scp[,c(3:50)])),]

pca_res_scp <- prcomp(t(pca_df_scp[,c(35:50)]))  # Transpose to have samples in rows
# Create a data frame with PCA results
pca_res_scp_df <- data.frame(
  PC1 = pca_res_scp$x[,1],
  PC2 = pca_res_scp$x[,2],
  condition = factor(c(rep("undiff",4),rep("diff",4),rep("aicar",4),rep("palm",4)))
)

# Calculate the proportion of variance explained
pct_var_scp_pca <- round(100 * (pca_res_scp$sdev^2 / sum(pca_res_scp$sdev^2)), digits = 2)

# Plot
pca_res_scp_df$label <- c(NA, "Undifferentiated\n+ Basal", rep(NA, 4), "Differentiated\n+ Basal",
                          rep(NA, 3), "Differentiated\n+ AICAR", rep(NA, 4), "Differentiated\n+ Palmitate")
pca_res_scp_df$color <- c(rep("#ff595e",4),rep("#ffca3a",4),rep("#8ac926",4), rep("#1982c4", 4))
pca_res_scp_df$status <- c(rep("undiff",4),rep("diff",12))
mpra_scp_pca <- ggplot(pca_res_scp_df, aes(PC1, PC2, color = condition, label=label, shape=status)) +
  geom_point(size=3, fill = pca_res_scp_df$color, color = "black") +
  scale_shape_manual(values = c(21, 24)) +
  xlab(paste0("PC1: ",pct_var_scp_pca[1],"% variance")) +
  ylab(paste0("PC2: ",pct_var_scp_pca[2],"% variance")) + 
  scale_color_manual(values = c("#ff595e","#ffca3a","#8ac926","#1982c4")) +
  geom_text_repel(min.segment.length = 0, family = "Helvetica",
                  box.padding = 1.5, size = 5, lineheight = 1,
                  data = pca_res_scp_df[c(1:4),], color = "black") + 
  geom_text_repel(min.segment.length = 0, family = "Helvetica",
                  box.padding = 1.5, nudge_y = 0.5, size = 5, lineheight = 1,
                  data = pca_res_scp_df[c(5:8),], color = "black") +
  geom_text_repel(min.segment.length = 0, family = "Helvetica",
                  box.padding = 1.5, nudge_x = -0.5, size = 5, lineheight = 1,
                  data = pca_res_scp_df[c(9:12),], color = "black") +
  geom_text_repel(min.segment.length = 0, family = "Helvetica",
                  box.padding = 1.5, nudge_y = -1, size = 5, lineheight = 1,
                  data = pca_res_scp_df[c(13:16),], color = "black") +
  theme_bw(base_family = "Helvetica", base_size = 16) + theme(aspect.ratio = 1, legend.position="none")

ggsave(mpra_scp_pca, filename = file.path(fig_dir, "mpra-scp-pca.png"),
       units = "in", dpi = 600, width = 5, height = 5, device = ragg::agg_png())

###############################################
### mpralm - Activity (Is an oligo active?) ###
###############################################

filt_cts <- filt_cts %>%
  mutate(category = case_when(str_detect(refname, "NA") ~ "neg_ctrl",
                              TRUE ~ "test_oligo")) %>%
  separate(refname, into = "rsid", remove = FALSE)

length(unique(filt_cts$rsid[filt_cts$prom=="MYBPC2"])) #294 vars tested
length(unique(filt_cts$rsid[filt_cts$prom=="SCP1"])) #290 vars tested

# Separate input DNA and RNA for mpralm formatting purposes
in_dna <- filt_cts %>%
  mutate(new_refname = paste(refname, prom, sep = "_")) %>%
  dplyr::select(starts_with("dna"), new_refname) %>%
  column_to_rownames(var = "new_refname")

in_rna <- filt_cts %>%
  mutate(new_refname = paste(refname, prom, sep = "_")) %>%
  dplyr::select(starts_with("rna"), new_refname) %>%
  column_to_rownames(var = "new_refname")

colnames(in_dna) <- colnames(in_rna) <- paste0("sample", seq(1,16,by=1))

# Perform activity analysis - basal undiff group
mpraset_undiff <- MPRASet(DNA = in_dna[,c(1:4)], RNA = in_rna[,c(1:4)], eid = row.names(in_dna))
design <- model.matrix(~1, data=data.frame(sample=paste0("sample", 1:4)))
undiff_fit <- mpralm(object=mpraset_undiff, design=design, aggregate="none",
                     model_type="indep_groups", plot=TRUE)

tr_undiff <- treat(undiff_fit)
active_res_undiff <- topTreat(tr_undiff, coef=1, number = Inf) %>%
  rownames_to_column(var = "new_refname")

write.table(active_res_undiff, file.path(out_dir, "active-res-undiff.tsv"),
            sep = "\t", quote = F, row.names = FALSE)
# diff group
mpraset_diff <- MPRASet(DNA = in_dna[,c(5:8)], RNA = in_rna[,c(5:8)], eid = row.names(in_dna))
diff_fit <- mpralm(object=mpraset_diff, design=design, aggregate="none",
                   model_type="indep_groups", plot=TRUE)

tr_diff <- treat(diff_fit)
active_res_diff <- topTreat(tr_diff, coef=1, number = Inf) %>%
  rownames_to_column(var = "new_refname")

write.table(active_res_diff, file.path(out_dir, "active-res-diff.tsv"),
            sep = "\t", quote = F, row.names = FALSE)
# AICAR group
mpraset_aicar <- MPRASet(DNA = in_dna[,c(9:12)], RNA = in_rna[,c(9:12)], eid = row.names(in_dna))
aicar_fit <- mpralm(object=mpraset_aicar, design=design, aggregate="none",
                    model_type="indep_groups", plot=TRUE)

tr_aicar <- treat(aicar_fit)
active_res_aicar <- topTreat(tr_aicar, coef=1, number = Inf) %>%
  rownames_to_column(var = "new_refname")

write.table(active_res_aicar, file.path(out_dir, "active-res-aicar.tsv"),
            sep = "\t", quote = F, row.names = FALSE)
# AICAR group
mpraset_palm <- MPRASet(DNA = in_dna[,c(13:16)], RNA = in_rna[,c(13:16)], eid = row.names(in_dna))
palm_fit <- mpralm(object=mpraset_palm, design=design, aggregate="none",
                   model_type="indep_groups", plot=TRUE)

tr_palm <- treat(palm_fit)
active_res_palm <- topTreat(tr_palm, coef=1, number = Inf) %>%
  rownames_to_column(var = "new_refname")

write.table(active_res_palm, file.path(out_dir, "active-res-palm.tsv"),
            sep = "\t", quote = F, row.names = FALSE)

# Filter to active (FDR < 0.05) oligos
lm_active_undiff = filtActive(active_res_undiff)
lm_active_diff = filtActive(active_res_diff)
lm_active_aicar = filtActive(active_res_aicar)
lm_active_palm = filtActive(active_res_palm)

# Filter to non-active (FDR > 0.05) oligos
lm_nonact_undiff = filtNonActive(active_res_undiff)
lm_nonact_diff = filtNonActive(active_res_diff)
lm_nonact_aicar = filtNonActive(active_res_aicar)
lm_nonact_palm = filtNonActive(active_res_palm)

# Filter to oligo allele pairs where at least one allele is active
lm_pairs_undiff <- processPairs(lm_active_undiff)
lm_pairs_diff <- processPairs(lm_active_diff)
lm_pairs_aicar <- processPairs(lm_active_aicar)
lm_pairs_palm <- processPairs(lm_active_palm)

lm_nonpairs_undiff <- processPairs(lm_nonact_undiff)
lm_nonpairs_diff <- processPairs(lm_nonact_diff)
lm_nonpairs_aicar <- processPairs(lm_nonact_aicar)
lm_nonpairs_palm <- processPairs(lm_nonact_palm)

lm_bck_undiff <- lm_nonpairs_undiff %>%
  filter(!(rsid %in% lm_pairs_undiff$rsid))
lm_bck_diff <- lm_nonpairs_diff %>%
  filter(!(rsid %in% lm_pairs_diff$rsid))
lm_bck_aicar <- lm_nonpairs_aicar %>%
  filter(!(rsid %in% lm_pairs_aicar$rsid))
lm_bck_palm <- lm_nonpairs_palm %>%
  filter(!(rsid %in% lm_pairs_palm$rsid))


##################################
### Make UpSet plot - Activity ###
##################################

# Make list of all active fragments from any group
lm_all_active_rsid = Reduce(union, list(lm_active_undiff$new_refname,
                                        lm_active_diff$new_refname,
                                        lm_active_aicar$new_refname,
                                        lm_active_palm$new_refname))

# Create individual boolean sets for each group
set_active_undiff = lm_all_active_rsid %in% lm_active_undiff$new_refname
set_active_diff = lm_all_active_rsid %in% lm_active_diff$new_refname
set_active_aicar = lm_all_active_rsid %in% lm_active_aicar$new_refname
set_active_palm = lm_all_active_rsid %in% lm_active_palm$new_refname

# Make final boolean input df
active_set = as.data.frame(cbind(set_active_undiff, set_active_diff,
                                 set_active_aicar, set_active_palm))
rownames(active_set) = lm_all_active_rsid

plot_names <- c("Undiff + Basal", "Diff + Basal",
                                        "Diff + AICAR", "Diff + Palmitate")

colnames(active_set) <- c("undiff", "diff", "aicar", "palm")

# Make list of all active rsids from any group
lm_all_active_pairs = Reduce(union, list(lm_pairs_undiff$new_refname,
                                        lm_pairs_diff$new_refname,
                                        lm_pairs_aicar$new_refname,
                                        lm_pairs_palm$new_refname))

# Create individual boolean sets for each group
pairs_active_undiff = lm_all_active_pairs %in% lm_pairs_undiff$new_refname
pairs_active_diff = lm_all_active_pairs %in% lm_pairs_diff$new_refname
pairs_active_aicar = lm_all_active_pairs %in% lm_pairs_aicar$new_refname
pairs_active_palm = lm_all_active_pairs %in% lm_pairs_palm$new_refname

# Make final boolean input df
active_pairs = as.data.frame(cbind(pairs_active_undiff, pairs_active_diff,
                                   pairs_active_aicar, pairs_active_palm))
rownames(active_pairs) = lm_all_active_pairs

colnames(active_pairs) <- c("undiff", "diff", "aicar", "palm")

active_set

# df for intersections
# - # fragments active in each set by promoter
active_intersect_df <- active_set %>%
  rownames_to_column(var = "new_refname") %>%
  separate(new_refname, into = c("refname", "promoter"), sep = "_(?=[^_]+$)") %>%
  filter(!str_detect(refname, "_NA$")) %>%
  group_by(refname, promoter) %>%
  arrange(promoter) %>%
  mutate(set_config = paste(undiff, diff, aicar, palm, sep = ","))

active_intersect_consist <- active_intersect_df %>%
  group_by(refname) %>%
  reframe(unique_prom = n_distinct(promoter), .groups = "drop",
          unique_configs = n_distinct(set_config), .groups = "drop") %>%
  mutate(is_consistent = unique_configs == 1)

active_intersect_df_consist <- active_intersect_df %>% left_join(active_intersect_consist, by = "refname")

active_intersect_plot_df <- active_intersect_df_consist %>%
  mutate(promoter_set = case_when(unique_prom == 2 & is_consistent == FALSE ~ promoter,
                                  unique_prom == 2 & is_consistent == TRUE ~ "both",
                                  unique_prom == 1 ~ promoter)) %>%
  ungroup() %>%
  select(-c(promoter, set_config, unique_configs, .groups, is_consistent))

active_upset_df <- upset(active_intersect_plot_df[,c(2:5,7)], 
                         intersect = c("undiff", "diff",
                                       "aicar", "palm"),
                         base_annotations=list(
                           'Intersection Size' = intersection_size(
                             counts = FALSE,
                             mapping = aes(linetype = promoter_set),
                             text = list(size = 5)
                           )
                         ),
                         plot_names,
                         queries = list(
                           upset_query(set = "palm", fill = "#1982c4"),
                           upset_query(set = "aicar", fill = "#8ac926"),
                           upset_query(set = "diff", fill = "#ffca3a"),
                           upset_query(set = "undiff", fill = "#ff595e")
                         ),
                         themes = upset_default_themes(
                           text = element_text(family = "Helvetica", size = 15)
                         ),
                         sort_sets = 'descending'
)

# df for all fragments per condition bars
# - # fragments active per condition, promoter
active_set_df <- active_set %>%
  rownames_to_column(var = "new_refname") %>%
  separate(new_refname, into = c("refname", "promoter"), sep = "_(?=[^_]+$)") %>%
  filter(!str_detect(refname, "_NA$")) %>%
  pivot_longer(!c(refname, promoter)) %>%
  group_by(refname, name) %>%
  reframe(refname = refname,
          promoter = paste(promoter, collapse = ","),
          value = paste(value, collapse = ",")) %>%
  filter(!(value %in% c("FALSE", "FALSE,FALSE"))) %>%
  ungroup() %>%
  distinct(.keep_all = TRUE) %>%
  mutate(promoter_set = case_when(
    promoter == "SCP1" & value == "TRUE" ~ "SCP1",
    promoter == "MYBPC2" & value == "TRUE" ~ "MYBPC2",
    promoter == "SCP1,MYBPC2" & value == "TRUE,TRUE" ~ "both",
    promoter == "SCP1,MYBPC2" & value == "TRUE,FALSE" ~ "SCP1",
    promoter == "SCP1,MYBPC2" & value == "FALSE,TRUE" ~ "MYBPC2",
    promoter == "MYBPC2,SCP1" & value == "TRUE,TRUE" ~ "both",
    promoter == "MYBPC2,SCP1" & value == "TRUE,FALSE" ~ "MYBPC2",
    promoter == "MYBPC2,SCP1" & value == "FALSE,TRUE" ~ "SCP1"))

# df for set size
active_set_tot <- active_set_df %>%
  group_by(promoter_set, name) %>%
  summarise(total = n_distinct(refname), .groups = "drop") %>%
  dplyr::rename(condition = name)

active_set_totals <- active_set_tot %>%
  group_by(condition) %>%
  reframe(label = sum(total))

cond_labels <- c("Undiff. +\nBasal","Diff. +\nBasal",
                 "Diff. +\nAICAR","Diff. +\nPalmitate"
)

active_set_bar <- active_set_tot %>% mutate(condition = factor(condition,
                                                               levels = c("undiff", "diff",
                                                                          "aicar", "palm"))) %>%
  ggplot(aes(x = condition, y = total, fill = condition, alpha = promoter_set)) +
  geom_bar(stat = "identity",
           linewidth = 0.25,
           color = "white",
           fill = "white",
           width = 0.6,
           alpha = 1) +
  geom_bar(stat = "identity",
           linewidth = 0.25,
           color = "black",
           width = 0.6) +
  geom_text(data = active_set_totals,
            aes(label = label, y = label + 5, x = condition),
            inherit.aes = FALSE,
            family = "Helvetica", vjust = 0, size = 4) +
  scale_alpha_manual(values = c(0.33, 0.66, 1),
                     labels = c("Both", "MYBPC2", "SCP1")) +
  scale_fill_manual(values = c("#ff595e","#ffca3a",
                               "#8ac926", "#1982c4")) +
  theme_linedraw(base_family = "Helvetica", base_size = 14) +
  theme(legend.position = c(0.8, 0.8),
        legend.spacing.y = unit(0.3, 'lines'),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.25),
        legend.key.size = unit(0.5,"line"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(color = "black", size = 14),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        plot.margin = margin(15, 0, 0, 0, "pt")) +
  labs(x = NULL, y = "Number of active oligos") + scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(labels = cond_labels) + coord_cartesian(clip = "off") +
  guides(alpha = guide_legend(byrow = TRUE, title = "Promoter", override.aes = list(label = "")),
         fill = "none")

ggsave(active_set_bar, filename = file.path(fig_dir, "active_set_bar.png"),
       units = "in", width = 5, height = 4, dpi = 600, device = ragg::agg_png())

# # of variants with at least one active allele
active_pairs_df <- active_pairs %>%
  rownames_to_column(var = "new_refname") %>%
  separate(new_refname, into = c("refname", "promoter"), sep = "_(?=[^_]+$)") %>%
  filter(!str_detect(refname, "neg_ctrl")) %>%
  pivot_longer(!c(refname, promoter)) %>%
  group_by(refname, name) %>%
  reframe(refname = refname,
          promoter = paste(promoter, collapse = ","),
          value = paste(value, collapse = ",")) %>%
  filter(!(value %in% c("FALSE", "FALSE,FALSE"))) %>%
  ungroup() %>%
  distinct(.keep_all = TRUE) %>%
  mutate(promoter_set = case_when(
    promoter == "SCP1" & value == "TRUE" ~ "SCP1",
    promoter == "MYBPC2" & value == "TRUE" ~ "MYBPC2",
    promoter == "SCP1,MYBPC2" & value == "TRUE,TRUE" ~ "both",
    promoter == "SCP1,MYBPC2" & value == "TRUE,FALSE" ~ "SCP1",
    promoter == "SCP1,MYBPC2" & value == "FALSE,TRUE" ~ "MYBPC2",
    promoter == "MYBPC2,SCP1" & value == "TRUE,TRUE" ~ "both",
    promoter == "MYBPC2,SCP1" & value == "TRUE,FALSE" ~ "MYBPC2",
    promoter == "MYBPC2,SCP1" & value == "FALSE,TRUE" ~ "SCP1"))

# df for set size
active_pairs_tot <- active_pairs_df %>%
  group_by(promoter_set, name) %>%
  summarise(total = n_distinct(refname), .groups = "drop") %>%
  dplyr::rename(condition = name)

active_pairs_totals <- active_pairs_tot %>%
  group_by(condition) %>%
  reframe(label = sum(total))

cond_labels <- c("Undiff. +\nBasal","Diff. +\nBasal",
                 "Diff. +\nAICAR","Diff. +\nPalmitate"
)


active_pairs_bar <- active_pairs_tot %>%
  mutate(condition = factor(condition,
                            levels = c("undiff", "diff",
                                       "aicar", "palm"))) %>%
  ggplot(aes(x = condition, y = total, fill = condition, alpha = promoter_set)) +
  geom_bar(stat = "identity",
           linewidth = 0.25,
           color = "#D3D3D3",
           fill = "#D3D3D3",
           width = 0.6,
           alpha = 1) +
  geom_bar(stat = "identity",
           linewidth = 0.25,
           color = "black",
           width = 0.6) +
  geom_text(data = active_pairs_totals,
            aes(label = label, y = label + 3, x = condition),
            inherit.aes = FALSE,
            family = "Helvetica", vjust = 0, size = 4) +
  scale_alpha_manual(values = c(0.25, 0.5, 0.75)) +
  scale_fill_manual(values = c("#ff595e","#ffca3a",
                               "#8ac926", "#1982c4")) +
  theme_linedraw(base_family = "Helvetica", base_size = 14) +
  theme(legend.position = "none",
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(color = "black", size = 14),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        plot.margin = margin(15, 0, 0, 0, "pt")) +
  labs(x = NULL, y = "Number of variants with\nat least 1 active allele") + scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(labels = cond_labels) + coord_cartesian(clip = "off")

ggsave(active_pairs_bar, filename = file.path(fig_dir, "active_pairs_bar.png"),
       units = "in", width = 5, height = 4, dpi = 600, device = ragg::agg_png())

# get levels of intersections
factor(active_upset_df[[2]]$data$intersection) # Levels: 1 1-3 2 2-1 2-1-3 2-3 3 4 4-1-3 4-2 4-2-1 4-2-1-3 4-2-3 4-3


active_intersect_df <- as.data.frame(table(factor(active_upset_df[[2]]$data$intersection),
                                           active_upset_df[[2]]$data$promoter_set))

active_intersect_totals <- active_intersect_df %>%
  group_by(Var1) %>%
  summarise(total = sum(Freq))
active_intersect_df <- active_intersect_df %>%
  left_join(active_intersect_totals, by = "Var1")
active_set_interbar <- active_intersect_df %>%
  mutate(condition = case_when(Var1 == 4 ~ "undiff",
                               Var1 == 2 ~ "diff",
                               Var1 == 1 ~ "aicar",
                               Var1 == 3 ~ "palm",
                               TRUE ~ "other"),
         condition = factor(condition, levels = c("undiff", "diff", "aicar", "palm", "other")),
         Var1 = reorder(Var1, -Freq)) %>%
  ggplot(aes(x = Var1, y = Freq, fill = condition, alpha = Var2)) +
  geom_bar(stat = "identity", linewidth = 0.25, color = "white", fill = "white", alpha = 1) +
  geom_bar(stat = "identity", linewidth = 0.25, color = "black") +
  geom_text(aes(label = total, y = total + 3),
            vjust = 0, size = 4, family = "Helvetica") +
  scale_alpha_manual(values = c(0.33, 0.66, 1)) + scale_fill_manual(values = c("#ff595e", "#ffca3a",
                                                                               "#8ac926", "#1982c4",
                                                                               "#8E8E8E")) +
  theme_minimal(base_family = "Helvetica", base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 14),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt")) +
  labs(x = NULL, y = "Intersection Size")

intersection_order <- levels(reorder(active_intersect_df$Var1, -active_intersect_df$Freq))
replace_intersect <- c("1" = "aicar", "2" = "diff", "3" = "palm", "4" = "undiff")
intersection_order <- as.factor(gsub("1", "aicar", 
                                     gsub("2", "diff", 
                                          gsub("3", "palm", 
                                               gsub("4", "undiff", intersection_order)))))
intersection_order <- factor(intersection_order, levels = intersection_order)

active_set_matrix <- active_upset_df[[4]]$data %>%
  mutate(start = str_sub(as.character(intersection), 1, 1),
         end = str_sub(as.character(intersection), nchar(as.character(intersection)), nchar(as.character(intersection))),
         intersection_recode = str_replace_all(intersection, replace_intersect),
         condition = case_when(intersection_recode == "undiff" ~ "undiff",
                               intersection_recode == "diff" ~ "diff",
                               intersection_recode == "aicar" ~ "aicar",
                               intersection_recode == "palm" ~ "palm",
                               TRUE ~ "other"),
         condition = factor(condition, levels = c("undiff", "diff", "aicar", "palm", "other")),
         group = str_replace_all(group, replace_intersect),
         group = factor(group, levels = c("palm", "aicar", "diff", "undiff")),
         start = str_replace_all(start, replace_intersect),
         end = str_replace_all(end, replace_intersect),
         intersection_recode = factor(intersection_recode, levels = intersection_order)) %>%
  ggplot(aes(x = intersection_recode, y = group, fill = condition,
             size = value)) +
  geom_point(fill = "#E4E4E2", color = "#C0C0C0", size = 4, pch = 21, stroke = 0.5) +
  geom_segment(aes(x = intersection_recode, xend = intersection_recode, y = start, yend = end),
               linewidth = 0.75) +
  geom_point(aes(stroke = value * 0.75), pch = 21) + scale_fill_manual(values = c("#ff595e", "#ffca3a",
                                                                                  "#8ac926", "#1982c4",
                                                                                  "#000000")) +
  scale_size_manual(values = c(0, 5)) + theme_minimal(base_family = "Helvetica", base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_text(color = "black", size = 14),
        plot.margin = margin(0, 0, 0, 0, "pt")) +
  labs(x = "Intersection", y = NULL) + scale_y_discrete(labels = c("Diff. + Palmitate","Diff. + AICAR",
                                                                   "Diff. + Basal", "Undiff. + Basal"))

design <- "
11
22
"
active_upset_plot <- free(active_set_interbar, type = "space", side = "l") + active_set_matrix + plot_layout(design = design)
active_bar_plots <- free(active_set_bar, type = "space", side = "l") + active_pairs_bar + plot_layout(design = design)

# save as pdf and add annotations to small intersections
ggsave(active_upset_plot, filename = file.path(fig_dir, "active_upset.pdf"),
       units = "in", width = 6, height = 6, dpi = 600)
ggsave(active_bar_plots, filename = file.path(fig_dir, "active_bars.png"),
       units = "in", width = 4, height = 6, dpi = 600)


########################
### mpralm - allelic ###
########################

# Add allele annotations and condense refnames
annot_filt_cts <- filt_cts %>%
  separate(refname, into = c("rsid", "chr", "pos", "ref", "alt",
                             "allele", "class", "site"),
           sep = "_", remove = FALSE) %>%
  mutate(class = case_when(ref == "NA" & alt == "NA" & allele == "NA" ~ "neg_ctrl",
                           alt == "plus" | alt == "minus" ~ "cage_strand",
                           class == "left" & is.na(site) ~ "standard_left",
                           class == "right" & is.na(site) ~ "standard_right",
                           site == "left" & class == "del" ~ "deletion_left",
                           site == "right" & class == "del" ~ "deletion_right",
                           site == "left" & class == "shf" ~ "shuffle_left",
                           site == "right" & class == "shf" ~ "shuffle_right",
                           class == "del" & is.na(site) ~ "deletion",
                           class == "shf" & is.na(site) ~ "shuffle",
                           TRUE ~ "standard"),
         allele = case_when(alt == "minus" ~ "A",
                            alt == "plus" ~ "R",
                            TRUE ~ allele),
         rsid = case_when(class == "cage_strand" ~ pos,
                          TRUE ~ rsid),
         rsid = paste(rsid, class, sep = "_"),
         new_refname = paste(rsid, prom, sep = "_"))

allele_filt_cts <- annot_filt_cts %>%
  dplyr::filter(!is.na(allele))

ref_cts <- allele_filt_cts %>%
  dplyr::filter(allele == "R")
alt_cts <- allele_filt_cts %>%
  dplyr::filter(allele == "A")

refnames <- colnames(ref_cts)[c(11:42)]
refnames <- paste(refnames, "ref", sep = "_")
colnames(ref_cts)[c(11:42)] <- refnames

altnames <- colnames(alt_cts)[c(11:42)]
altnames <- paste(altnames, "alt", sep = "_")
colnames(alt_cts)[c(11:42)] <- altnames

allelic_cts <- merge(ref_cts[,c(2:6,10:42,44)], alt_cts[,c(2:6,10:42,44)],
                     by=c('rsid'='rsid','chr'='chr',
                          'pos'='pos','ref'='ref',
                          'alt'='alt','prom'='prom',
                          'new_refname'='new_refname'))

allelic_info <- allelic_cts[,c(1:7)]
allelic_info <- allelic_info %>%
  extract(new_refname, into = c("refname", "prom"), "(.*)_([^_]+)$")

write.table(allelic_info, 
            file.path(out_dir, "allelic-info.tsv"),
            sep = "\t", quote = F, row.names = F)

# Now filter data to include only oligo pairs where at least one oligo is active
# Start with basal group
in_dna_allelic <- allelic_cts %>%
  dplyr::select(new_refname, starts_with("dna"))
in_rna_allelic <- allelic_cts %>%
  dplyr::select(new_refname, starts_with("rna"))

# prep for mpralm
in_dna_undiff <- in_dna_allelic %>%
  dplyr::select(dna1_ref:dna4_ref, dna1_alt:dna4_alt, new_refname) %>%
  dplyr::filter(new_refname %in% lm_pairs_undiff$new_refname) %>%
  column_to_rownames(var = "new_refname")

in_rna_undiff <- in_rna_allelic %>%
  dplyr::select(rna1_ref:rna4_ref,rna1_alt:rna4_alt, new_refname) %>%
  dplyr::filter(new_refname %in% lm_pairs_undiff$new_refname) %>%
  column_to_rownames(var = "new_refname")

colnames(in_dna_undiff) <- colnames(in_rna_undiff) <-
  c(paste(paste0("sample", seq(1,4,by=1)), "ref", sep = "_"),
    paste(paste0("sample", seq(1,4,by=1)), "alt", sep = "_"))

# mpralm allelic testing
mpraset_undiff_allele <- MPRASet(DNA = in_dna_undiff, RNA = in_rna_undiff,
                                 eid = row.names(in_dna_undiff))
# set up design matrix; each column is one of 4 treatment groups, each row is a sample
design = data.frame(intcpt = 1,
                    "alleleALT" = c(rep(0,4),rep(1,4)))
row.names(design) <- colnames(in_dna_undiff)
block_vector <- rep(1:4, 2)

undiff_fit_allele = mpralm(mpraset_undiff_allele, design=design, aggregate="none",
                           model_type="corr_groups",
                           block=block_vector, plot=TRUE)
# create summary table
allele_res_undiff <- topTable(undiff_fit_allele, coef = 2, number = Inf) %>%
  rownames_to_column(var = "new_refname")

write.table(allele_res_undiff, file.path(out_dir, "allele-res-undiff.tsv"),
            sep = "\t", quote = F, row.names = F)

# prep for mpralm - diff
in_dna_diff <- in_dna_allelic %>%
  dplyr::select(dna5_ref:dna8_ref, dna5_alt:dna8_alt, new_refname) %>%
  dplyr::filter(new_refname %in% lm_pairs_diff$new_refname) %>%
  column_to_rownames(var = "new_refname")

in_rna_diff <- in_rna_allelic %>%
  dplyr::select(rna5_ref:rna8_ref,rna5_alt:rna8_alt, new_refname) %>%
  dplyr::filter(new_refname %in% lm_pairs_diff$new_refname) %>%
  column_to_rownames(var = "new_refname")

colnames(in_dna_diff) <- colnames(in_rna_diff) <-
  c(paste(paste0("sample", seq(1,4,by=1)), "ref", sep = "_"),
    paste(paste0("sample", seq(1,4,by=1)), "alt", sep = "_"))

# mpralm allelic testing
mpraset_diff_allele <- MPRASet(DNA = in_dna_diff, RNA = in_rna_diff,
                               eid = row.names(in_dna_diff))
# set up design matrix; each column is one of 4 treatment groups, each row is a sample
design = data.frame(intcpt = 1,
                    "alleleALT" = c(rep(0,4),rep(1,4)))
row.names(design) <- colnames(in_dna_diff)
block_vector <- rep(1:4, 2)

diff_fit_allele = mpralm(mpraset_diff_allele, design=design, aggregate="none",
                         model_type="corr_groups",
                         block=block_vector, plot=TRUE)
# create summary table
allele_res_diff <- topTable(diff_fit_allele, coef = 2, number = Inf) %>%
  rownames_to_column(var = "new_refname")

write.table(allele_res_diff, file.path(out_dir, "allele-res-diff.tsv"),
            sep = "\t", quote = F, row.names = F)

# prep for mpralm - aicar
in_dna_aicar <- in_dna_allelic %>%
  dplyr::select(dna9_ref:dna12_ref, dna9_alt:dna12_alt, new_refname) %>%
  dplyr::filter(new_refname %in% lm_pairs_aicar$new_refname) %>%
  column_to_rownames(var = "new_refname")

in_rna_aicar <- in_rna_allelic %>%
  dplyr::select(rna9_ref:rna12_ref,rna9_alt:rna12_alt, new_refname) %>%
  dplyr::filter(new_refname %in% lm_pairs_aicar$new_refname) %>%
  column_to_rownames(var = "new_refname")

colnames(in_dna_aicar) <- colnames(in_rna_aicar) <-
  c(paste(paste0("sample", seq(1,4,by=1)), "ref", sep = "_"),
    paste(paste0("sample", seq(1,4,by=1)), "alt", sep = "_"))

# mpralm allelic testing
mpraset_aicar_allele <- MPRASet(DNA = in_dna_aicar, RNA = in_rna_aicar,
                                eid = row.names(in_dna_aicar))
# set up design matrix; each column is one of 4 treatment groups, each row is a sample
design = data.frame(intcpt = 1,
                    "alleleALT" = c(rep(0,4),rep(1,4)))
row.names(design) <- colnames(in_dna_aicar)
block_vector <- rep(1:4, 2)

aicar_fit_allele = mpralm(mpraset_aicar_allele, design=design, aggregate="none",
                          model_type="corr_groups",
                          block=block_vector, plot=TRUE)
# create summary table
allele_res_aicar <- topTable(aicar_fit_allele, coef = 2, number = Inf) %>%
  rownames_to_column(var = "new_refname")

write.table(allele_res_aicar, file.path(out_dir, "allele-res-aicar.tsv"),
            sep = "\t", quote = F, row.names = F)

# prep for mpralm - palm
in_dna_palm <- in_dna_allelic %>%
  dplyr::select(dna13_ref:dna16_ref, dna13_alt:dna16_alt, new_refname) %>%
  dplyr::filter(new_refname %in% lm_pairs_palm$new_refname) %>%
  column_to_rownames(var = "new_refname")

in_rna_palm <- in_rna_allelic %>%
  dplyr::select(rna13_ref:rna16_ref,rna13_alt:rna16_alt, new_refname) %>%
  dplyr::filter(new_refname %in% lm_pairs_palm$new_refname) %>%
  column_to_rownames(var = "new_refname")

colnames(in_dna_palm) <- colnames(in_rna_palm) <-
  c(paste(paste0("sample", seq(1,4,by=1)), "ref", sep = "_"),
    paste(paste0("sample", seq(1,4,by=1)), "alt", sep = "_"))

# mpralm allelic testing
mpraset_palm_allele <- MPRASet(DNA = in_dna_palm, RNA = in_rna_palm,
                               eid = row.names(in_dna_palm))
# set up design matrix; each column is one of 4 treatment groups, each row is a sample
design = data.frame(intcpt = 1,
                    "alleleALT" = c(rep(0,4),rep(1,4)))
row.names(design) <- colnames(in_dna_palm)
block_vector <- rep(1:4, 2)

palm_fit_allele = mpralm(mpraset_palm_allele, design=design, aggregate="none",
                         model_type="corr_groups",
                         block=block_vector, plot=TRUE)
# create summary table
allele_res_palm <- topTable(palm_fit_allele, coef = 2, number = Inf) %>%
  rownames_to_column(var = "new_refname")

write.table(allele_res_palm, file.path(out_dir, "allele-res-palm.tsv"),
            sep = "\t", quote = F, row.names = F)

#################################
### Make UpSet plot - Allelic ###
#################################

# Get allelic results
lm_allele_undiff <- allele_res_undiff %>%
  extract(new_refname, into = c("refname", "prom"), "(.*)_([^_]+)$") %>%
  dplyr::filter(adj.P.Val < 0.1)

lm_nonall_undiff <- allele_res_undiff %>%
  extract(new_refname, into = c("refname", "prom"), "(.*)_([^_]+)$") %>%
  dplyr::filter(adj.P.Val > 0.1)

lm_allele_diff <- allele_res_diff %>%
  extract(new_refname, into = c("refname", "prom"), "(.*)_([^_]+)$") %>%
  dplyr::filter(adj.P.Val < 0.1)

lm_nonall_diff <- allele_res_diff %>%
  extract(new_refname, into = c("refname", "prom"), "(.*)_([^_]+)$") %>%
  dplyr::filter(adj.P.Val > 0.1)

lm_allele_aicar <- allele_res_aicar %>%
  extract(new_refname, into = c("refname", "prom"), "(.*)_([^_]+)$") %>%
  dplyr::filter(adj.P.Val < 0.1)

lm_nonall_aicar <- allele_res_aicar %>%
  extract(new_refname, into = c("refname", "prom"), "(.*)_([^_]+)$") %>%
  dplyr::filter(adj.P.Val > 0.1)

lm_allele_palm <- allele_res_palm %>%
  extract(new_refname, into = c("refname", "prom"), "(.*)_([^_]+)$") %>%
  dplyr::filter(adj.P.Val < 0.1)

lm_nonall_palm <- allele_res_palm %>%
  extract(new_refname, into = c("refname", "prom"), "(.*)_([^_]+)$") %>%
  dplyr::filter(adj.P.Val > 0.1)


rsid_allele_undiff <- paste(lm_allele_undiff$refname, lm_allele_undiff$prom, sep = "_")
rsid_allele_diff <- paste(lm_allele_diff$refname, lm_allele_diff$prom, sep = "_")
rsid_allele_aicar <- paste(lm_allele_aicar$refname, lm_allele_aicar$prom, sep = "_")
rsid_allele_palm <- paste(lm_allele_palm$refname, lm_allele_palm$prom, sep = "_")

rsid_allele_set = Reduce(union, list(rsid_allele_undiff,rsid_allele_diff,
                                     rsid_allele_aicar,rsid_allele_palm))
undiff_rsid_set = rsid_allele_set %in% rsid_allele_undiff
diff_rsid_set = rsid_allele_set %in% rsid_allele_diff
aicar_rsid_set = rsid_allele_set %in% rsid_allele_aicar
palm_rsid_set = rsid_allele_set %in% rsid_allele_palm

allelic_set = as.data.frame(cbind(undiff_rsid_set, diff_rsid_set,
                                  aicar_rsid_set,palm_rsid_set))
rownames(allelic_set) <- rsid_allele_set
colnames(allelic_set) <- c("undiff","diff","aicar","palm")

# also make set w logfc
rsid_lfc_undiff <- lm_allele_undiff %>%
  mutate(new_refname = paste(refname, prom, sep = "_")) %>%
  filter(new_refname %in% rsid_allele_set) %>%
  select(new_refname, logFC)
rsid_lfc_diff <- lm_allele_diff %>%
  mutate(new_refname = paste(refname, prom, sep = "_")) %>%
  filter(new_refname %in% rsid_allele_set) %>%
  select(new_refname, logFC)
rsid_lfc_aicar <- lm_allele_aicar %>%
  mutate(new_refname = paste(refname, prom, sep = "_")) %>%
  filter(new_refname %in% rsid_allele_set) %>%
  select(new_refname, logFC)
rsid_lfc_palm <- lm_allele_palm %>%
  mutate(new_refname = paste(refname, prom, sep = "_")) %>%
  filter(new_refname %in% rsid_allele_set) %>%
  select(new_refname, logFC)

lfc_set <- allelic_set %>%
  rownames_to_column(var = "new_refname") %>%
  full_join(rsid_lfc_undiff) %>%
  rename(undiff_lfc = logFC) %>%
  full_join(rsid_lfc_diff) %>%
  rename(diff_lfc = logFC) %>%
  full_join(rsid_lfc_aicar) %>%
  rename(aicar_lfc = logFC) %>%
  full_join(rsid_lfc_palm) %>%
  rename(palm_lfc = logFC) %>%
  rowwise() %>%
  mutate(consistent = case_when(
    all(c_across(undiff_lfc:palm_lfc) < 0, na.rm = TRUE) ~ TRUE,
    all(c_across(undiff_lfc:palm_lfc) > 0, na.rm = TRUE) ~ TRUE,
    TRUE ~ FALSE)
  ) %>%
  ungroup()

binom_res <- lfc_set %>%
  summarise(
    consistent_count = sum(consistent, na.rm = TRUE),
    total_count = n()
  ) %>%
  mutate(
    binomial_test = list(rstatix::binom_test(consistent_count, total_count, p = 0.5))
  ) %>%
  unnest_wider(binomial_test)


write.table(lfc_set, file.path(out_dir, "lfc_allele_comparison.tsv"),
                               sep = "\t", quote = F, row.names = F)

allelic_set_df <- allelic_set %>%
  rownames_to_column(var = "new_refname") %>%
  separate(new_refname, into = c("refname", "promoter"), sep = "_(?=[^_]+$)") %>%
  pivot_longer(!c(refname, promoter)) %>%
  group_by(refname, name) %>%
  reframe(refname = refname,
          promoter = paste(promoter, collapse = ","),
          value = paste(value, collapse = ",")) %>%
  filter(!(value %in% c("FALSE", "FALSE,FALSE"))) %>%
  ungroup() %>%
  distinct(.keep_all = TRUE) %>%
  mutate(promoter_set = case_when(
    promoter == "SCP1" & value == "TRUE" ~ "SCP1",
    promoter == "MYBPC2" & value == "TRUE" ~ "MYBPC2",
    promoter == "SCP1,MYBPC2" & value == "TRUE,TRUE" ~ "both",
    promoter == "SCP1,MYBPC2" & value == "TRUE,FALSE" ~ "SCP1",
    promoter == "SCP1,MYBPC2" & value == "FALSE,TRUE" ~ "MYBPC2",
    promoter == "MYBPC2,SCP1" & value == "TRUE,TRUE" ~ "both",
    promoter == "MYBPC2,SCP1" & value == "TRUE,FALSE" ~ "MYBPC2",
    promoter == "MYBPC2,SCP1" & value == "FALSE,TRUE" ~ "SCP1"))

allelic_intersect_df <- allelic_set %>%
  rownames_to_column(var = "new_refname") %>%
  separate(new_refname, into = c("refname", "promoter"), sep = "_(?=[^_]+$)") %>%
  filter(!str_detect(refname, "_NA$")) %>%
  group_by(refname, promoter) %>%
  arrange(promoter) %>%
  mutate(set_config = paste(undiff, diff, aicar, palm, sep = ","))

allelic_intersect_consist <- allelic_intersect_df %>%
  group_by(refname) %>%
  reframe(unique_prom = n_distinct(promoter), .groups = "drop",
          unique_configs = n_distinct(set_config), .groups = "drop") %>%
  mutate(is_consistent = unique_configs == 1)

allelic_intersect_df_consist <- allelic_intersect_df %>% left_join(allelic_intersect_consist, by = "refname")

allelic_intersect_plot_df <- allelic_intersect_df_consist %>%
  mutate(promoter_set = case_when(unique_prom == 2 & is_consistent == FALSE ~ promoter,
                                  unique_prom == 2 & is_consistent == TRUE ~ "both",
                                  unique_prom == 1 ~ promoter)) %>%
  ungroup() %>%
  select(-c(promoter, set_config, unique_configs, .groups, is_consistent))


allelic_upset_df <- upset(allelic_intersect_plot_df[,c(2:5,7)], 
                         intersect = c("undiff", "diff",
                                       "aicar", "palm"),
                         base_annotations=list(
                           'Intersection Size' = intersection_size(
                             counts = FALSE,
                             mapping = aes(linetype = promoter_set),
                             text = list(size = 5)
                           )
                         ),
                         plot_names,
                         queries = list(
                           upset_query(set = "palm", fill = "#1982c4"),
                           upset_query(set = "aicar", fill = "#8ac926"),
                           upset_query(set = "diff", fill = "#ffca3a"),
                           upset_query(set = "undiff", fill = "#ff595e")
                         ),
                         themes = upset_default_themes(
                           text = element_text(family = "Helvetica", size = 15)
                         ),
                         sort_sets = 'descending'
)

allelic_set_tot <- allelic_set_df %>%
  group_by(promoter_set, name) %>%
  summarise(total = n_distinct(refname), .groups = "drop") %>%
  dplyr::rename(condition = name)

allelic_set_totals <- allelic_set_tot %>%
  group_by(condition) %>%
  summarise(label = sum(total))
allelic_set_bar <- allelic_set_tot %>% mutate(condition = factor(condition,
                                                                     levels = c("undiff", "diff",
                                                                                "aicar", "palm"))) %>%
  ggplot(aes(x = condition, y = total, fill = condition, alpha = promoter_set)) +
  geom_bar(stat = "identity",
           linewidth = 0.25,
           color = "white",
           fill = "white",
           width = 0.6,
           alpha = 1) +
  geom_bar(stat = "identity",
           linewidth = 0.25,
           color = "black",
           width = 0.6) +
  geom_text(data = allelic_set_totals,
            aes(label = label, y = label + 4, x = condition),
            inherit.aes = FALSE,
            family = "Helvetica", vjust = 0, size = 4) +
  scale_alpha_manual(values = c(0.33, 0.66, 1),
                     labels = c("Both", "MYBPC2", "SCP1")) +
  scale_fill_manual(values = c("#ff595e","#ffca3a",
                               "#8ac926", "#1982c4")) +
  theme_linedraw(base_family = "Helvetica", base_size = 14) +
  theme(legend.position = c(0.8, 0.85),
        legend.spacing.y = unit(0.3, 'lines'),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.25),
        legend.key.size = unit(0.5,"line"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(color = "black", size = 14),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        plot.margin = margin(20, 0, 0, 0, "pt")) +
  labs(x = NULL, y = "Number of variants with allelic bias") + scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(labels = cond_labels) + coord_cartesian(clip = "off") +
  guides(alpha = guide_legend(byrow = TRUE, title = "Promoter", override.aes = list(label = "")),
         fill = "none")

ggsave(allelic_set_bar, filename = file.path(fig_dir, "allelic_set_bar.png"),
       units = "in", width = 5, height = 4, dpi = 600, device = ragg::agg_png())


# get levels of intersections
factor(allelic_upset_df[[2]]$data$intersection) # Levels: 1 1-3 2 2-1 2-1-3 2-3 3 4 4-1 4-1-3 4-2 4-2-1 4-2-1-3 4-2-3 4-3


allelic_intersect_df <- as.data.frame(table(factor(allelic_upset_df[[2]]$data$intersection),
                                            allelic_upset_df[[2]]$data$promoter_set))

allelic_intersect_totals <- allelic_intersect_df %>%
  group_by(Var1) %>%
  summarise(total = sum(Freq))
allelic_intersect_df <- allelic_intersect_df %>%
  left_join(allelic_intersect_totals, by = "Var1")
allelic_set_interbar <- allelic_intersect_df %>%
  mutate(condition = case_when(Var1 == 4 ~ "undiff",
                               Var1 == 2 ~ "diff",
                               Var1 == 1 ~ "aicar",
                               Var1 == 3 ~ "palm",
                               TRUE ~ "other"),
         condition = factor(condition, levels = c("undiff", "diff", "aicar", "palm", "other")),
         Var1 = reorder(Var1, -Freq)) %>%
  ggplot(aes(x = Var1, y = Freq, fill = condition, alpha = Var2)) +
  geom_bar(stat = "identity", linewidth = 0.25, color = "white", fill = "white", alpha = 1) +
  geom_bar(stat = "identity", linewidth = 0.25, color = "black") +
  geom_text(aes(label = total, y = total + 3),
            vjust = 0, size = 4, family = "Helvetica") +
  scale_alpha_manual(values = c(0.33, 0.66, 1),
                     labels = c("Both", "MYBPC2", "SCP1")) + scale_fill_manual(values = c("#ff595e", "#ffca3a",
                                                                                          "#8ac926", "#1982c4",
                                                                                          "#8E8E8E")) +
  theme_minimal(base_family = "Helvetica", base_size = 14) +
  theme(legend.position = c(0.85, 0.75),
        legend.spacing.y = unit(0.3, 'lines'),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.25),
        legend.key.size = unit(0.5,"line"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 14),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt")) +
  labs(x = NULL, y = "Intersection Size") +
  guides(alpha = guide_legend(byrow = TRUE, title = "Promoter", override.aes = list(label = "")),
         fill = "none")

intersection_order <- levels(reorder(allelic_intersect_df$Var1, -allelic_intersect_df$Freq))
replace_intersect <- c("1" = "aicar", "2" = "diff", "3" = "palm", "4" = "undiff")
intersection_order <- as.factor(gsub("1", "aicar", 
                                     gsub("2", "diff", 
                                          gsub("3", "palm", 
                                               gsub("4", "undiff", intersection_order)))))
intersection_order <- factor(intersection_order, levels = intersection_order)

allelic_set_matrix <- allelic_upset_df[[4]]$data %>%
  mutate(start = str_sub(as.character(intersection), 1, 1),
         end = str_sub(as.character(intersection), nchar(as.character(intersection)), nchar(as.character(intersection))),
         intersection_recode = str_replace_all(intersection, replace_intersect),
         condition = case_when(intersection_recode == "undiff" ~ "undiff",
                               intersection_recode == "diff" ~ "diff",
                               intersection_recode == "aicar" ~ "aicar",
                               intersection_recode == "palm" ~ "palm",
                               TRUE ~ "other"),
         condition = factor(condition, levels = c("undiff", "diff", "aicar", "palm", "other")),
         group = str_replace_all(group, replace_intersect),
         group = factor(group, levels = c("palm", "aicar", "diff", "undiff")),
         start = str_replace_all(start, replace_intersect),
         end = str_replace_all(end, replace_intersect),
         intersection_recode = factor(intersection_recode, levels = intersection_order)) %>%
  ggplot(aes(x = intersection_recode, y = group, fill = condition,
             size = value)) +
  geom_point(fill = "#E4E4E2", color = "#C0C0C0", size = 4, pch = 21, stroke = 0.5) +
  geom_segment(aes(x = intersection_recode, xend = intersection_recode, y = start, yend = end),
               linewidth = 0.75) +
  geom_point(aes(stroke = value * 0.75), pch = 21) + scale_fill_manual(values = c("#ff595e", "#ffca3a",
                                                                                  "#8ac926", "#1982c4",
                                                                                  "#000000")) +
  scale_size_manual(values = c(0, 5)) + theme_minimal(base_family = "Helvetica", base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_text(color = "black", size = 14),
        plot.margin = margin(0, 0, 0, 0, "pt")) +
  labs(x = "Intersection", y = NULL) + scale_y_discrete(labels = c("Diff. + Palmitate","Diff. + AICAR",
                                                                   "Diff. + Basal", "Undiff. + Basal"))

design <- "
11
11
11
22
22
"

allelic_upset_plot <- free(allelic_set_interbar, type = "space", side = "l") + allelic_set_matrix + plot_layout(design = design)

# save as pdf and add annotations to small intersections
ggsave(allelic_upset_plot, filename = file.path(fig_dir, "allelic_upset.pdf"),
       units = "in", width = 6, height = 6, dpi = 600)


############################
### Single var plot fxns ###
############################

# Function to manipulate dfs for p-val testing and plotting
make_rsid_df <- function(neg_ids, ref_ids, plot_prom, df) {
  
  # make neg control df
  neg_df <- df %>% dplyr::select(-c(rna1:dna16)) %>%
    dplyr::filter(refname %in% neg_ids & prom == plot_prom) %>%
    dplyr::select(-prom) %>%
    pivot_longer(!c(refname),
                 names_to = "replicate",
                 values_to = "value") %>%
    group_by(replicate) %>%
    reframe(value = mean(value), .groups = 'drop') %>%
    mutate(refname = "neg_ctrl")
  
  mean_neg_value <- neg_df %>%
    summarize(mean_value = mean(value)) %>%
    pull(mean_value)
  
  test_df <- df %>% dplyr::select(-c(rna1:dna16)) %>%
    dplyr::filter(refname %in% ref_ids & prom == plot_prom) %>%
    dplyr::select(-prom) %>%
    pivot_longer(!c(refname),
                 names_to = "replicate",
                 values_to = "value") %>%
    bind_rows(neg_df) %>%
    dplyr::filter(is.finite(value)) %>%
    mutate(plot_value = value - mean_neg_value) %>%
    mutate(allele = case_when(refname == "neg_ctrl" ~ "NC",
                              substr(refname, nchar(refname)-1, nchar(refname)) == "ft" ~ substr(refname, nchar(refname)-5,nchar(refname)-5),
                              substr(refname, nchar(refname)-1, nchar(refname)) == "ht" ~ substr(refname, nchar(refname)-6,nchar(refname)-6),
                              TRUE ~ substr(refname, nchar(refname), nchar(refname))),
           allele = factor(allele, levels = c("NC", "R", "A")),
           condition = case_when(replicate %in% c("ratio1","ratio2",
                                                  "ratio3","ratio4") ~ "undiff",
                                 replicate %in% c("ratio5", "ratio6",
                                                  "ratio7","ratio8") ~ "diff",
                                 replicate %in% c("ratio9","ratio10",
                                                  "ratio11","ratio12") ~ "aicar",
                                 replicate %in% c("ratio13","ratio14",
                                                  "ratio15","ratio16") ~ "palm"),
           condition = factor(condition, levels = c("undiff", "diff", "aicar","palm")))
  
  if(length(test_df$replicate) %% 4 == 0){
    pval_df <- test_df %>% group_by(condition) %>%
      rstatix::t_test(value ~ allele, paired = TRUE) %>%
      rstatix::add_xy_position(x = "allele", dodge = 0.75) %>%
      mutate(p.format = case_when(p < 0.005 ~ paste0("p = ", scientific(p, digits = 3)),
                                  TRUE ~ paste0("p = ", p)))
  } else {
    pval_df <- test_df %>% group_by(condition) %>%
      rstatix::t_test(value ~ allele) %>%
      rstatix::add_xy_position(x = "allele", dodge = 0.75) %>%
      mutate(p.format = case_when(p < 0.005 ~ paste0("p = ", scientific(p, digits = 3)),
                                  TRUE ~ paste0("p = ", p)))
  }
  
  result <- setNames(
    list(neg_df, test_df, pval_df),
    c("neg_df", "test_df", "pval_df")
  )
  
  return(result)
  
}

# Function to make rsid facet plot separated by condition
make_rsid_plot <- function(test_df, pval_df, rsid, allele_labels, condition_colors, condition_labeller) {
  
  test_df %>% ggplot(aes(x = allele, y = plot_value)) +
    geom_hline(yintercept = 0, linewidth = 0.15) +
    geom_boxplot(aes(color = factor(condition)), outlier.shape = NA,
                 position = position_dodge()) + 
    geom_point(aes(fill = factor(condition)),
               color = "black", pch = 21, size = 2, stroke = 0.5,
               position = position_dodge(width = 0.75)) +
    facet_wrap(~condition, strip.position="bottom",
               labeller = labeller(condition = condition_labeller),
               nrow = 1) +
    labs(x = NULL, y = bquote(log[2]~"(RNA/DNA), " ~ .(rsid))) +
    scale_x_discrete(labels = allele_labels) +
    theme_bw(base_family = "Helvetica", base_size = 16) +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_blank(),
          axis.ticks = element_blank(),
          axis.title.y = element_text(size = 14),
          strip.placement = "outside",
          strip.background = element_blank(),
          strip.text = element_text(size = 16)) +
    scale_fill_manual(values = condition_colors) +
    scale_color_manual(values = condition_colors) +
    stat_pvalue_manual(
      pval_df, label = "p.format", tip.length = 0.01,
      family = "Helvetica", size = 4.5)
  
}

#################
### Plot vars ###
#################

condition_colors <- c("#ff595e","#ffca3a",
                      "#8ac926","#1982c4")
condition_labeller <- c(
  undiff = "Undiff. +\nBasal",
  diff = "Diff. +\nBasal",
  aicar = "Diff. +\nAICAR",
  palm = "Diff. +\nPalmitate"
)

rs118_out <- make_rsid_df(neg_ids = c("NA_chr5_116670618_NA_NA_NA"),
                          ref_ids = c("rs11867290_chr17_62791047_T_G_A",
                                      "rs11867290_chr17_62791047_T_G_R"),
                          plot_prom = "MYBPC2", df = plot_tpm)

rs118_out$pval_df <- rs118_out$pval_df %>%
  mutate(y.position = c(0.6, 0.75, 1, 1.15, 1.3, 1.45,
                        1.15, 1.3, 1.45, 1.15, 1.3, 1.45))

rs118_plot <- make_rsid_plot(test_df = rs118_out$test_df,
                             pval_df = rs118_out$pval_df,
                             rsid = "rs11867290",
                             allele_labels = c("Neg. Ctrl","Ref (T)", "Alt (G)"),
                             condition_colors, condition_labeller) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  scale_y_continuous(limits = c(-0.25, 1.5)) + labs(y = bquote(log[2]~"(RNA/DNA), " ~ .("rs11867290-MYBPC2")))

rs118_plot
ggsave(rs118_plot, filename = file.path(fig_dir, "rs11867290-plot.png"),
       units = "in", dpi = 600, width = 7.2, height = 5.5)

rs340_out <- make_rsid_df(neg_ids = c("NA_chr5_116670618_NA_NA_NA"),
                          ref_ids = c("rs34003091_chr19_58150887_T_C_R",
                                      "rs34003091_chr19_58150887_T_C_A"),
                          plot_prom = "MYBPC2", df = plot_tpm)

rs340_out$pval_df <- rs340_out$pval_df %>%
  mutate(y.position = c(1.65, 1.85, 2.05, 1, 1.2, 1.4,
                        0.9, 1.1, 1.3, 0.9, 1.1, 1.3))

rs340_plot <- make_rsid_plot(test_df = rs340_out$test_df,
                             pval_df = rs340_out$pval_df,
                             rsid = "rs34003091",
                             allele_labels = c("Neg. Ctrl","Ref (T)", "Alt (C)"),
                             condition_colors, condition_labeller) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  scale_y_continuous(limits = c(-0.25, 2.1)) + labs(y = bquote(log[2]~"(RNA/DNA), " ~ .("rs34003091-MYBPC2")))

rs340_plot
ggsave(rs340_plot, filename = file.path(fig_dir, "rs34003091-plot.png"),
       units = "in", dpi = 600, width = 7.2, height = 5.5)

rs490_out <- make_rsid_df(neg_ids = "NA_chr15_74506364_NA_NA_NA",
                          ref_ids = c("rs490972_chr11_66312315_G_A_A",
                                      "rs490972_chr11_66312315_G_A_R"),
                          plot_prom = "SCP1", df = plot_tpm)

rs490_out$pval_df <- rs490_out$pval_df %>%
  mutate(y.position = c(0.825, 0.95, 1.075, 0.25, 0.375, 0.5,
                        0.2, 0.325, 0.475, 0.1, 0.225, 0.375))

rs490_plot <- make_rsid_plot(test_df = rs490_out$test_df,
                             pval_df = rs490_out$pval_df,
                             rsid = "rs490972",
                             allele_labels = c("Neg. Ctrl","Ref (G)", "Alt (A)"),
                             condition_colors, condition_labeller) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  scale_y_continuous(limits = c(-0.55, 1.1)) + labs(y = bquote(log[2]~"(RNA/DNA), " ~ .("rs490972-SCP1")))


rs490_plot
ggsave(rs490_plot, filename = file.path(fig_dir, "rs490972-plot.png"),
       units = "in", dpi = 600, width = 9.6, height = 5)
write.table(rs490_out$test_df, file.path(out_dir, "rs490972_test-df.tsv"), sep = "\t", quote = F, row.names = F)


#################
### ISGU data ###
#################

formatted_isgu <- read.csv("/lab/work/knishino/20240731_lhcn-isgu/20240731_formatted-lhcn-isgu-aggregate-names-updated.csv")
background_isgu <- read.csv("/lab/work/knishino/20240731_lhcn-isgu/20240731_formatted-lhcn-isgu - background_wells_not_adjacent.csv")


formatted_isgu <- formatted_isgu %>%
  rename(sample = Sample, condition = Condition, insulin = Insulin,
         thirty = X30min, sixty = X60min, ninety = X90min, onetwenty = X120min) %>%
  pivot_longer(!c(sample, condition, insulin)) 

background_isgu <- background_isgu %>%
  rename(thirty = X30min, sixty = X60min, ninety = X90min, onetwenty = X120min) %>%
  pivot_longer(thirty:onetwenty) %>%
  group_by(name) %>%
  summarize(avg_background = mean(value))

normalized_isgu <- formatted_isgu %>%
  left_join(background_isgu, by = "name") %>%
  filter(!(sample %in% c("ATS0666","ATS0673"))) %>%
  mutate(norm_value = value - avg_background,
         condition = case_when(condition == "basal" ~ "diff",
                               condition == "AICAR" ~ "aicar",
                               condition == "palmitate" ~ "palm",
                               TRUE ~ condition),
         condition = factor(condition, levels = c("undiff", "diff",
                                                  "aicar","palm")))

pval_df_ins <- normalized_isgu %>% filter(name == "sixty") %>% group_by(condition) %>%
  rstatix::t_test(norm_value ~ insulin) %>%
  rstatix::add_xy_position(x = "condition", dodge = 0.75) %>%
  mutate(p.format = case_when(p < 0.005 ~ paste0("p = ", scientific(p, digits = 3)),
                              TRUE ~ paste0("p = ", p)),
         y.position = c(2.2e6,
                        8.1e6,
                        7.8e6,
                        7e6))

pval_df_cond <- normalized_isgu %>% filter(name == "sixty" & insulin == FALSE) %>%
  rstatix::t_test(norm_value ~ condition) %>%
  rstatix::add_xy_position(x = "condition", dodge = 0.75) %>%
  mutate(p.format = case_when(p < 0.005 ~ paste0("p = ", scientific(p, digits = 3)),
                              TRUE ~ paste0("p = ", p)))
pval_df_cond <- pval_df_cond[c(1:3),] %>%
  mutate(y.position = c(2.75e6,3.25e6,3.75e6))

norm_isgu_plot <- normalized_isgu %>%
  filter(name == "sixty") %>%
  ggplot(aes(x = condition, y = norm_value)) +
  geom_boxplot(aes(fill = condition, group = interaction(insulin, condition)), outlier.shape = NA,
               position = position_dodge()) +
  geom_point(aes(shape = insulin),
             color = "black", fill = "white", size = 3, stroke = 0.5,
             position = position_jitterdodge()) +
  labs(x = NULL, y = "Relative Light Units (RLU)") +
  scale_x_discrete(labels = condition_labeller) +
  theme_bw(base_family = "Helvetica", base_size = 16) +
  theme(legend.position = c(0.15, 0.85),
        legend.spacing.y = unit(0.3, 'lines'),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.25),
        legend.key.size = unit(0.5,"line"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14)) +
  scale_fill_manual(values = condition_colors) +
  scale_shape_manual(values = c(21,19),
                     labels = c("- insulin",
                                "+ insulin"),
                     name = "Stimulation") +
  scale_y_continuous(limits = c(0, 8.5e6),
                     labels = c(0, expression(paste("2 x ", 10^6)),
                                expression(paste("4 x ", 10^6)),
                                expression(paste("6 x ", 10^6)),
                                expression(paste("8 x ", 10^6)))) +
  stat_pvalue_manual(
    pval_df_ins, label = "p.format", tip.length = 0.01,
    family = "Helvetica", size = 4.5
  ) +
  stat_pvalue_manual(
    pval_df_cond, label = "p.format", tip.length = 0.01,
    family = "Helvetica", size = 4.5
  ) +
  guides(fill="none")

ggsave(norm_isgu_plot, filename = file.path(fig_dir, "norm-isgu-plot.png"),
       units = "in", dpi = 600, width = 7, height = 5)
