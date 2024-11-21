#############################
### Install/load packages ###
#############################

#install.packages("librarian")
librarian::shelf(tidyverse, factoextra, data.table, mpra, reshape2,
                 ggtext, UpSetR, patchwork, GGally, ComplexUpset,
                 EnsDb.Hsapiens.v79, locuszoomr, ensembldb, Biostrings,
                 ggrepel, gkmSVM, rstatix, ggpubr, scales, BSgenome.Hsapiens.UCSC.hg38,
                 chromVARmotifs, ggh4x)

#####################################
### Functions for processing data ###
#####################################

filtActive <- function(res){
  res %>% dplyr::filter(adj.P.Val < 0.05) %>% dplyr::select(refname)
}

filtNonActive <- function(res){
  res %>% dplyr::filter(adj.P.Val > 0.05) %>% dplyr::select(refname)
}

processPairs = function(df) {
  df %>%
    extract(refname, into = c("refname", "prom"), "(.*)_([^_]+)$") %>%
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
  geom_text_repel(min.segment.length = 0, family = "Helvetica",
                  box.padding = 1.5, size = 5, lineheight = 1,
                  data = pca_res_df[c(1:4),], color = "black") + 
  geom_text_repel(min.segment.length = 0, family = "Helvetica",
                  box.padding = 1.5, nudge_y = -0.05, size = 5, lineheight = 1,
                  data = pca_res_df[c(5:8),], color = "black") +
  geom_text_repel(min.segment.length = 0, family = "Helvetica",
                  box.padding = 1.5, nudge_y = -0.1, nudge_x = -0.05, size = 5, lineheight = 1,
                  data = pca_res_df[c(9:12),], color = "black") +
  geom_text_repel(min.segment.length = 0, family = "Helvetica",
                  box.padding = 1.5, nudge_x = -0.05, size = 5, lineheight = 1,
                  data = pca_res_df[c(13:16),], color = "black") +
  theme_bw(base_family = "Helvetica", base_size = 16) + theme(aspect.ratio = 1, legend.position="none")

ggsave(mpra_pca, filename = file.path(fig_dir, "mpra-pca.png"),
       units = "in", dpi = 600, width = 5, height = 5, device = ragg::agg_png())

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
                  box.padding = 1.5, nudge_y = -0.1, nudge_x = -0.05, size = 5, lineheight = 1,
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
                  box.padding = 1.5, nudge_y = -0.05, size = 5, lineheight = 1,
                  data = pca_res_scp_df[c(5:8),], color = "black") +
  geom_text_repel(min.segment.length = 0, family = "Helvetica",
                  box.padding = 1.5, nudge_y = -0.35, nudge_x = -0.05, size = 5, lineheight = 1,
                  data = pca_res_scp_df[c(9:12),], color = "black") +
  geom_text_repel(min.segment.length = 0, family = "Helvetica",
                  box.padding = 1.5, nudge_x = -0.05, size = 5, lineheight = 1,
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
            sep = "\t", quote = F)
# diff group
mpraset_diff <- MPRASet(DNA = in_dna[,c(5:8)], RNA = in_rna[,c(5:8)], eid = row.names(in_dna))
diff_fit <- mpralm(object=mpraset_diff, design=design, aggregate="none",
                   model_type="indep_groups", plot=TRUE)

tr_diff <- treat(diff_fit)
active_res_diff <- topTreat(tr_diff, coef=1, number = Inf) %>%
  rownames_to_column(var = "new_refname")

write.table(active_res_diff, file.path(out_dir, "active-res-diff.tsv"),
            sep = "\t", quote = F)
# AICAR group
mpraset_aicar <- MPRASet(DNA = in_dna[,c(9:12)], RNA = in_rna[,c(9:12)], eid = row.names(in_dna))
aicar_fit <- mpralm(object=mpraset_aicar, design=design, aggregate="none",
                    model_type="indep_groups", plot=TRUE)

tr_aicar <- treat(aicar_fit)
active_res_aicar <- topTreat(tr_aicar, coef=1, number = Inf) %>%
  rownames_to_column(var = "new_refname")

write.table(active_res_aicar, file.path(out_dir, "active-res-aicar.tsv"),
            sep = "\t", quote = F)
# AICAR group
mpraset_palm <- MPRASet(DNA = in_dna[,c(13:16)], RNA = in_rna[,c(13:16)], eid = row.names(in_dna))
palm_fit <- mpralm(object=mpraset_palm, design=design, aggregate="none",
                   model_type="indep_groups", plot=TRUE)

tr_palm <- treat(palm_fit)
active_res_palm <- topTreat(tr_palm, coef=1, number = Inf) %>%
  rownames_to_column(var = "new_refname")

write.table(active_res_palm, file.path(out_dir, "active-res-palm.tsv"),
            sep = "\t", quote = F)

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

##################################
### Make UpSet plot - Activity ###
##################################

# Make list of all active rsids from any group
lm_all_active_rsid = Reduce(union, list(lm_pairs_undiff$new_refname,
                                        lm_pairs_diff$new_refname,
                                        lm_pairs_aicar$new_refname,
                                        lm_pairs_palm$new_refname))

# Create individual boolean sets for each group
set_active_undiff = lm_all_active_rsid %in% lm_pairs_undiff$new_refname
set_active_diff = lm_all_active_rsid %in% lm_pairs_diff$new_refname
set_active_aicar = lm_all_active_rsid %in% lm_pairs_aicar$new_refname
set_active_palm = lm_all_active_rsid %in% lm_pairs_palm$new_refname

# Make final boolean input df
active_set = as.data.frame(cbind(set_active_undiff, set_active_diff,
                                 set_active_aicar, set_active_palm))
rownames(active_set) = lm_all_active_rsid

plot_names <- colnames(active_set) <- c("Undiff + Basal", "Diff + Basal",
                                        "Diff + AICAR", "Diff + Palmitate")

png(file.path(fig_dir, "upset-active.png"), width = 8, height = 4.5, units = "in",
    res = 300)
upset(active_set, plot_names, name = "Intersection",
      queries = list(
        upset_query(intersect = "Diff + Palmitate", color = "#1982c4", fill = "#1982c4"),
        upset_query(intersect = "Diff + AICAR", color = "#8ac926", fill = "#8ac926"),
        upset_query(intersect = "Diff + Basal", color = "#ffca3a", fill = "#ffca3a"),
        upset_query(intersect = "Undiff + Basal", color = "#ff595e", fill = "#ff595e"),
        upset_query(set = "Diff + Palmitate", fill = "#1982c4"),
        upset_query(set = "Diff + AICAR", fill = "#8ac926"),
        upset_query(set = "Diff + Basal", fill = "#ffca3a"),
        upset_query(set = "Undiff + Basal", fill = "#ff595e")
      ),
      themes = upset_default_themes(
        text = element_text(family = "Helvetica", size = 15)
      ),
      sort_sets = 'descending',
      base_annotations = list(
        'Intersection size' = intersection_size(
          text = list(size = 5)
        )
      )
)
dev.off()

###############################################
### mpralm - Activity (Is an oligo active?)
### MYBPC2
###############################################

filt_cts_myb <- filt_cts %>%
  dplyr::filter(prom == "MYBPC2")

# Separate input DNA and RNA for mpralm formatting purposes
in_dna <- filt_cts_myb %>%
  mutate(new_refname = paste(refname, prom, sep = "_")) %>%
  dplyr::select(starts_with("dna"), new_refname) %>%
  column_to_rownames(var = "new_refname")

in_rna <- filt_cts_myb %>%
  mutate(new_refname = paste(refname, prom, sep = "_")) %>%
  dplyr::select(starts_with("rna"), new_refname) %>%
  column_to_rownames(var = "new_refname")

colnames(in_dna) <- colnames(in_rna) <- paste0("sample", seq(1,16,by=1))

# Perform activity analysis - basal undiff group
mpraset_undiff <- MPRASet(DNA = in_dna[,c(1:4)], RNA = in_rna[,c(1:4)], eid = row.names(in_dna))
design <- model.matrix(~1, data=data.frame(sample=paste0("sample", 1:4)))
undiff_fit_myb <- mpralm(object=mpraset_undiff, design=design, aggregate="none",
                         model_type="indep_groups", plot=TRUE)

tr_undiff_myb <- treat(undiff_fit_myb)
active_res_undiff_myb <- topTreat(tr_undiff_myb, coef=1, number = Inf) %>%
  rownames_to_column(var = "new_refname")

write.table(active_res_undiff_myb, file.path(out_dir, "active-res-undiff-myb.tsv"),
            sep = "\t", quote = F)
# diff group
mpraset_diff <- MPRASet(DNA = in_dna[,c(5:8)], RNA = in_rna[,c(5:8)], eid = row.names(in_dna))
diff_fit_myb <- mpralm(object=mpraset_diff, design=design, aggregate="none",
                       model_type="indep_groups", plot=TRUE)

tr_diff_myb <- treat(diff_fit)
active_res_diff_myb <- topTreat(tr_diff_myb, coef=1, number = Inf) %>%
  rownames_to_column(var = "new_refname")

write.table(active_res_diff_myb, file.path(out_dir, "active-res-diff-myb.tsv"),
            sep = "\t", quote = F)
# AICAR group
mpraset_aicar <- MPRASet(DNA = in_dna[,c(9:12)], RNA = in_rna[,c(9:12)], eid = row.names(in_dna))
aicar_fit_myb <- mpralm(object=mpraset_aicar, design=design, aggregate="none",
                        model_type="indep_groups", plot=TRUE)

tr_aicar_myb <- treat(aicar_fit_myb)
active_res_aicar_myb <- topTreat(tr_aicar_myb, coef=1, number = Inf) %>%
  rownames_to_column(var = "new_refname")

write.table(active_res_aicar_myb, file.path(out_dir, "active-res-aicar-myb.tsv"),
            sep = "\t", quote = F)
# AICAR group
mpraset_palm <- MPRASet(DNA = in_dna[,c(13:16)], RNA = in_rna[,c(13:16)], eid = row.names(in_dna))
palm_fit_myb <- mpralm(object=mpraset_palm, design=design, aggregate="none",
                       model_type="indep_groups", plot=TRUE)

tr_palm_myb <- treat(palm_fit_myb)
active_res_palm_myb <- topTreat(tr_palm_myb, coef=1, number = Inf) %>%
  rownames_to_column(var = "new_refname")

write.table(active_res_palm_myb, file.path(out_dir, "active-res-palm-myb.tsv"),
            sep = "\t", quote = F)

# Filter to active (FDR < 0.05) oligos
lm_active_undiff_myb = filtActive(active_res_undiff_myb)
lm_active_diff_myb = filtActive(active_res_diff_myb)
lm_active_aicar_myb = filtActive(active_res_aicar_myb)
lm_active_palm_myb = filtActive(active_res_palm_myb)

# Filter to oligo allele pairs where at least one allele is active
lm_pairs_undiff_myb <- processPairs(lm_active_undiff_myb)
lm_pairs_diff_myb <- processPairs(lm_active_diff_myb)
lm_pairs_aicar_myb <- processPairs(lm_active_aicar_myb)
lm_pairs_palm_myb <- processPairs(lm_active_palm_myb)

##################################
### Make UpSet plot - Activity ###
##################################

# Make list of all active rsids from any group
lm_all_active_rsid_myb = Reduce(union, list(lm_pairs_undiff_myb$new_refname,
                                            lm_pairs_diff_myb$new_refname,
                                            lm_pairs_aicar_myb$new_refname,
                                            lm_pairs_palm_myb$new_refname))

# Create individual boolean sets for each group
set_active_undiff_myb = lm_all_active_rsid_myb %in% lm_pairs_undiff_myb$new_refname
set_active_diff_myb = lm_all_active_rsid_myb %in% lm_pairs_diff_myb$new_refname
set_active_aicar_myb = lm_all_active_rsid_myb %in% lm_pairs_aicar_myb$new_refname
set_active_palm_myb = lm_all_active_rsid_myb %in% lm_pairs_palm_myb$new_refname

# Make final boolean input df
active_set_myb = as.data.frame(cbind(set_active_undiff_myb, set_active_diff_myb,
                                     set_active_aicar_myb, set_active_palm_myb))


plot_names <- colnames(active_set_myb) <- c("Undiff + Basal", "Diff + Basal",
                                            "Diff + AICAR", "Diff + Palmitate")

png(file.path(fig_dir, "upset-active-myb.png"), width = 8, height = 4.5, units = "in",
    res = 300)
upset(active_set_myb, plot_names, name = "Intersection",
      queries = list(
        upset_query(intersect = "Diff + Palmitate", color = "#1982c4", fill = "#1982c4"),
        upset_query(intersect = "Diff + AICAR", color = "#8ac926", fill = "#8ac926"),
        upset_query(intersect = "Diff + Basal", color = "#ffca3a", fill = "#ffca3a"),
        upset_query(intersect = "Undiff + Basal", color = "#ff595e", fill = "#ff595e"),
        upset_query(set = "Diff + Palmitate", fill = "#1982c4"),
        upset_query(set = "Diff + AICAR", fill = "#8ac926"),
        upset_query(set = "Diff + Basal", fill = "#ffca3a"),
        upset_query(set = "Undiff + Basal", fill = "#ff595e")
      ),
      themes = upset_default_themes(
        text = element_text(family = "Helvetica", size = 15)
      ),
      sort_sets = 'descending',
      base_annotations = list(
        'Intersection size' = intersection_size(
          text = list(size = 5)
        )
      )
)
dev.off()

###############################################
### mpralm - Activity (Is an oligo active?) ###
###############################################

filt_cts_scp <- filt_cts %>%
  dplyr::filter(prom == "SCP1")

# Separate input DNA and RNA for mpralm formatting purposes
in_dna <- filt_cts_scp %>%
  mutate(new_refname = paste(refname, prom, sep = "_")) %>%
  dplyr::select(starts_with("dna"), new_refname) %>%
  column_to_rownames(var = "new_refname")

in_rna <- filt_cts_scp %>%
  mutate(new_refname = paste(refname, prom, sep = "_")) %>%
  dplyr::select(starts_with("rna"), new_refname) %>%
  column_to_rownames(var = "new_refname")

colnames(in_dna) <- colnames(in_rna) <- paste0("sample", seq(1,16,by=1))

# Perform activity analysis - basal undiff group
mpraset_undiff <- MPRASet(DNA = in_dna[,c(1:4)], RNA = in_rna[,c(1:4)], eid = row.names(in_dna))
design <- model.matrix(~1, data=data.frame(sample=paste0("sample", 1:4)))
undiff_fit_scp <- mpralm(object=mpraset_undiff, design=design, aggregate="none",
                         model_type="indep_groups", plot=TRUE)

tr_undiff_scp <- treat(undiff_fit_scp)
active_res_undiff_scp <- topTreat(tr_undiff_scp, coef=1, number = Inf) %>%
  rownames_to_column(var = "new_refname")

write.table(active_res_undiff_scp, file.path(out_dir, "active-res-undiff-scp.tsv"),
            sep = "\t", quote = F)
# diff group
mpraset_diff <- MPRASet(DNA = in_dna[,c(5:8)], RNA = in_rna[,c(5:8)], eid = row.names(in_dna))
diff_fit_scp <- mpralm(object=mpraset_diff, design=design, aggregate="none",
                       model_type="indep_groups", plot=TRUE)

tr_diff_scp <- treat(diff_fit_scp)
active_res_diff_scp <- topTreat(tr_diff_scp, coef=1, number = Inf) %>%
  rownames_to_column(var = "new_refname")

write.table(active_res_diff_scp, file.path(out_dir, "active-res-diff-scp.tsv"),
            sep = "\t", quote = F)
# AICAR group
mpraset_aicar <- MPRASet(DNA = in_dna[,c(9:12)], RNA = in_rna[,c(9:12)], eid = row.names(in_dna))
aicar_fit_scp <- mpralm(object=mpraset_aicar, design=design, aggregate="none",
                        model_type="indep_groups", plot=TRUE)

tr_aicar_scp <- treat(aicar_fit_scp)
active_res_aicar_scp <- topTreat(tr_aicar_scp, coef=1, number = Inf) %>%
  rownames_to_column(var = "new_refname")

write.table(active_res_aicar_scp, file.path(out_dir, "active-res-aicar_scp.tsv"),
            sep = "\t", quote = F)
# AICAR group
mpraset_palm <- MPRASet(DNA = in_dna[,c(13:16)], RNA = in_rna[,c(13:16)], eid = row.names(in_dna))
palm_fit_scp <- mpralm(object=mpraset_palm, design=design, aggregate="none",
                       model_type="indep_groups", plot=TRUE)

tr_palm_scp <- treat(palm_fit_scp)
active_res_palm_scp <- topTreat(tr_palm_scp, coef=1, number = Inf) %>%
  rownames_to_column(var = "new_refname")

write.table(active_res_palm_scp, file.path(out_dir, "active-res-palm_scp.tsv"),
            sep = "\t", quote = F)

# Filter to active (FDR < 0.05) oligos
lm_active_undiff_scp = filtActive(active_res_undiff_scp)
lm_active_diff_scp = filtActive(active_res_diff_scp)
lm_active_aicar_scp = filtActive(active_res_aicar_scp)
lm_active_palm_scp = filtActive(active_res_palm_scp)

# Filter to oligo allele pairs where at least one allele is active
lm_pairs_undiff_scp <- processPairs(lm_active_undiff_scp)
lm_pairs_diff_scp <- processPairs(lm_active_diff_scp)
lm_pairs_aicar_scp <- processPairs(lm_active_aicar_scp)
lm_pairs_palm_scp <- processPairs(lm_active_palm_scp)

##################################
### Make UpSet plot - Activity ###
##################################

# Make list of all active rsids from any group
lm_all_active_rsid_scp = Reduce(union, list(lm_pairs_undiff_scp$new_refname,
                                            lm_pairs_diff_scp$new_refname,
                                            lm_pairs_aicar_scp$new_refname,
                                            lm_pairs_palm_scp$new_refname))

# Create individual boolean sets for each group
set_active_undiff_scp = lm_all_active_rsid_scp %in% lm_pairs_undiff_scp$new_refname
set_active_diff_scp = lm_all_active_rsid_scp %in% lm_pairs_diff_scp$new_refname
set_active_aicar_scp = lm_all_active_rsid_scp %in% lm_pairs_aicar_scp$new_refname
set_active_palm_scp = lm_all_active_rsid_scp %in% lm_pairs_palm_scp$new_refname

# Make final boolean input df
active_set_scp = as.data.frame(cbind(set_active_undiff_scp, set_active_diff_scp,
                                     set_active_aicar_scp, set_active_palm_scp))


plot_names <- colnames(active_set_scp) <- c("Undiff + Basal", "Diff + Basal",
                                            "Diff + AICAR", "Diff + Palmitate")

png(file.path(fig_dir, "upset-active-scp.png"), width = 8, height = 4.5, units = "in",
    res = 300)
upset(active_set_scp, plot_names, name = "Intersection",
      themes = upset_default_themes(
        text = element_text(family = "Helvetica", size = 15)
      ),
      sort_sets = 'descending',
      base_annotations = list(
        'Intersection size' = intersection_size(
          text = list(size = 5)
        )
      )
)
dev.off()


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
            sep = "\t", quote = F)

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
            sep = "\t", quote = F)

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
            sep = "\t", quote = F)

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
            sep = "\t", quote = F)

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

colnames(allelic_set) <- plot_names

png(file.path(fig_dir, "upset-allelic.png"),
    width = 8, height = 4.5, units = "in",res = 300)
upset(allelic_set, plot_names, name = "Intersection",
      queries = list(
        upset_query(intersect = "Diff + Palmitate", color = "#1982c4", fill = "#1982c4"),
        upset_query(intersect = "Diff + AICAR", color = "#8ac926", fill = "#8ac926"),
        upset_query(intersect = "Diff + Basal", color = "#ffca3a", fill = "#ffca3a"),
        upset_query(intersect = "Undiff + Basal", color = "#ff595e", fill = "#ff595e"),
        upset_query(set = "Diff + Palmitate", fill = "#1982c4"),
        upset_query(set = "Diff + AICAR", fill = "#8ac926"),
        upset_query(set = "Diff + Basal", fill = "#ffca3a"),
        upset_query(set = "Undiff + Basal", fill = "#ff595e")
      ),
      themes = upset_default_themes(
        text = element_text(family = "Helvetica", size = 15)
      ),
      sort_sets = 'descending',
      base_annotations = list(
        'Intersection size' = intersection_size(
          text = list(size = 5)
        )
      )
)
dev.off()

############################
### Single var plot fxns ###
############################

#rs34003091

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

# options: rs181704186_chr1_159752293_A_G_R
# rs11419307_chr7_99520455_TTTTTTTTTT_TTTTTTTTTTT_R

rs114_out <- make_rsid_df(neg_ids = c("NA_chr11_116619500_NA_NA_NA",
                                      "NA_chr10_72885496_NA_NA_NA"),
                          ref_ids = c("rs11419307_chr7_99520455_TTTTTTTTTT_TTTTTTTTTTT_R",
                                      "rs11419307_chr7_99520455_TTTTTTTTTT_TTTTTTTTTTT_A"),
                          plot_prom = "MYBPC2", df = plot_tpm)

rs114_out$pval_df <- rs114_out$pval_df %>%
  mutate(y.position = c(2, 2.25, 2.5, 2.25, 2.5, 2.75,
                        1.75, 2, 2.25, 1.75, 2, 2.25))

condition_colors <- c("#ff595e","#ffca3a",
                      "#8ac926","#1982c4")
condition_labeller <- c(
  undiff = "Undiff. +\nBasal",
  diff = "Diff. +\nBasal",
  aicar = "Diff. +\nAICAR",
  palm = "Diff. +\nPalmitate"
)

rs114_plot <- make_rsid_plot(test_df = rs114_out$test_df,
                             pval_df = rs114_out$pval_df,
                             rsid = "rs11419307",
                             allele_labels = c("Neg. Ctrl","Ref (T)10", "Alt (T)11"),
                             condition_colors, condition_labeller) +
  scale_y_continuous(limits = c(-0.85, 2.85))

rs114_plot
ggsave(rs114_plot, filename = file.path(fig_dir, "rs11419307-plot.png"),
       units = "in", dpi = 600, width = 10, height = 4.5)


rs181_out <- make_rsid_df(neg_ids = c("NA_chr11_116619500_NA_NA_NA"),
                          ref_ids = c("rs181704186_chr1_159752293_A_G_R",
                                      "rs181704186_chr1_159752293_A_G_A"),
                          plot_prom = "SCP1", df = plot_tpm)

rs181_out$pval_df <- rs181_out$pval_df %>%
  mutate(y.position = c(1.15, 1.4, 1.65, 0.25, 0.5, 0.75,
                        0.25, 0.5, 0.75, 0.25, 0.5, 0.75))

rs181_plot <- make_rsid_plot(test_df = rs181_out$test_df,
                             pval_df = rs181_out$pval_df,
                             rsid = "rs181704186",
                             allele_labels = c("Neg. Ctrl","Ref (A)", "Alt (G)"),
                             condition_colors, condition_labeller) +
  scale_y_continuous(limits = c(-0.9, 1.75))

rs181_plot
ggsave(rs181_plot, filename = file.path(fig_dir, "rs181704186-plot.png"),
       units = "in", dpi = 600, width = 10, height = 4.5)

rs381_out <- make_rsid_df(neg_ids = "NA_chr10_72885496_NA_NA_NA",
                          ref_ids = c("rs3810155_chr19_10654191_C_G_R",
                                      "rs3810155_chr19_10654191_C_G_A"),
                          plot_prom = "SCP1", df = plot_tpm)

rs381_out$pval_df <- rs381_out$pval_df %>%
  mutate(y.position = c(2.25, 2.5, 2.75, 1.2, 1.45, 1.7,
                        0.4, 0.65, 0.9, 0.85, 1.1, 1.35))

rs381_plot <- make_rsid_plot(test_df = rs381_out$test_df,
                             pval_df = rs381_out$pval_df,
                             rsid = "rs3810155",
                             allele_labels = c("Neg. Ctrl","Ref (C)", "Alt (G)"),
                             condition_colors, condition_labeller) +
  scale_y_continuous(limits = c(-0.25, 2.8))

rs381_plot
ggsave(rs381_plot, filename = file.path(fig_dir, "rs3810155-plot.png"),
       units = "in", dpi = 600, width = 10, height = 4.5)

rs490_out <- make_rsid_df(neg_ids = "NA_chr10_72885496_NA_NA_NA",
                          ref_ids = c("rs490972_chr11_66312315_G_A_A",
                                      "rs490972_chr11_66312315_G_A_R"),
                          plot_prom = "SCP1", df = plot_tpm)

rs490_out$pval_df <- rs490_out$pval_df %>%
  mutate(y.position = c(0.825, 0.95, 1.075, 0.25, 0.375, 0.5,
                        0.25, 0.375, 0.5, 0.19, 0.315, 0.44))

rs490_plot <- make_rsid_plot(test_df = rs490_out$test_df,
                             pval_df = rs490_out$pval_df,
                             rsid = "rs490972",
                             allele_labels = c("Neg. Ctrl","Ref (G)", "Alt (A)"),
                             condition_colors, condition_labeller) +
  scale_y_continuous(limits = c(-0.55, 1.1))

rs490_plot
ggsave(rs490_plot, filename = file.path(fig_dir, "rs490972-plot.png"),
       units = "in", dpi = 600, width = 10, height = 4.5)
write.table(rs490_out$test_df, file.path(out_dir, "rs490972_test-df.tsv"), sep = "\t", quote = F, row.names = F)

rs313_out <- make_rsid_df(neg_ids = "NA_chr10_72885496_NA_NA_NA",
                          ref_ids = c("rs3130288_chr6_32128224_C_A_A",
                                      "rs3130288_chr6_32128224_C_A_R"),
                          plot_prom = "MYBPC2", df = plot_tpm)

rs381_out$pval_df <- rs381_out$pval_df %>%
  mutate(y.position = c(2.25, 2.5, 2.75, 1.2, 1.45, 1.7,
                        0.4, 0.65, 0.9, 0.85, 1.1, 1.35))

rs313_plot <- make_rsid_plot(test_df = rs313_out$test_df,
                             pval_df = rs313_out$pval_df,
                             rsid = "rs3130288",
                             allele_labels = c("Neg. Ctrl","Ref (C)", "Alt (A)"),
                             condition_colors, condition_labeller)

rs313_plot
ggsave(rs381_plot, filename = file.path(fig_dir, "rs3810155-plot.png"),
       units = "in", dpi = 600, width = 10, height = 4.5)
write.table(rs381_out$test_df, file.path(out_dir, "rs3810155_test-df.tsv"), sep = "\t", quote = F, row.names = F)

rs468_out <- make_rsid_df(neg_ids = "NA_chr10_72885496_NA_NA_NA",
                          ref_ids = c("rs4689391_chr4_6278722_G_A_A",
                                      "rs4689391_chr4_6278722_G_A_R"),
                          plot_prom = "MYBPC2", df = plot_tpm)

rs468_out$pval_df <- rs468_out$pval_df %>%
  mutate(y.position = c(2.25, 2.5, 2.75, 1.2, 1.45, 1.7,
                        0.4, 0.65, 0.9, 0.85, 1.1, 1.35))

rs468_plot <- make_rsid_plot(test_df = rs468_out$test_df,
                             pval_df = rs468_out$pval_df,
                             rsid = "rs4689391",
                             allele_labels = c("Neg. Ctrl","Ref (G)", "Alt (A)"),
                             condition_colors, condition_labeller)

+
  scale_y_continuous(limits = c(-0.25, 2.8))

rs381_plot
ggsave(rs381_plot, filename = "20241119_rs3810155-plot.png",
       units = "in", dpi = 600, width = 10, height = 4.5)


###################################
### compare effects across prom ###
###################################

compare_basal <- active_res_undiff_myb %>%
  rownames_to_column(var = "refname") %>%
  extract(refname, into = c("refname", "prom"), "(.*)_([^_]+)$") %>%
  inner_join(
    active_res_undiff_scp %>%
      rownames_to_column(var = "refname") %>%
      extract(refname, into = c("refname", "prom"), "(.*)_([^_]+)$"),
    by = "refname"
  ) %>%
  mutate(signed_pval_myb = case_when(logFC.x < 0 ~ -adj.P.Val.x,
                                     TRUE ~ adj.P.Val.x),
         signed_pval_scp = case_when(logFC.y < 0 ~ - adj.P.Val.y,
                                     TRUE ~ adj.P.Val.y))
compare_basal %>%
  filter(abs(signed_pval_scp) < 0.5 & abs(signed_pval_myb) < 0.5) %>%
  ggplot(aes(x = signed_pval_scp, y = signed_pval_myb)) +
  geom_point()

########################
### Motif enrichment ###
########################

options(meme_bin = "/lab/sw/modules/meme/5.0.4/bin")

# Activity sets 
get_active_seq <- function(df, df_null, genome) {
  
  active_coord <- df %>%
    separate(refname, into = c("rsid", "chr", "pos")) %>%
    distinct(rsid, .keep_all = TRUE) %>%
    dplyr::filter(rsid != "NA") %>%
    mutate(start = as.numeric(pos) - 100,
           stop = as.numeric(pos) + 99) %>%
    distinct(chr, start, stop, .keep_all = TRUE)
  
  active_coord_gr <- makeGRangesFromDataFrame(active_coord)
  
  nonact_coord <- df_null %>%
    separate(refname, into = c("rsid", "chr", "pos")) %>%
    distinct(rsid, .keep_all = TRUE) %>%
    dplyr::filter(rsid != "NA") %>%
    mutate(start = as.numeric(pos) - 100,
           stop = as.numeric(pos) + 99) %>%
    distinct(chr, start, stop, .keep_all = TRUE)
  
  nonact_coord_gr <- makeGRangesFromDataFrame(nonact_coord)
  
  active_seq <- active_coord_gr %>%
    get_sequence(genome)
  
  nonactive_seq = nonact_coord_gr %>%
    get_sequence(genome)
  
  return(list(active_seq = active_seq,
              nonactive_seq = nonactive_seq))
  
}

hg.genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
jaspar_path <- file.path(in_dir, "jaspar.meme")

undiff_active_seq <- get_active_seq(lm_active_undiff, lm_nonact_undiff, hg.genome)
undiff_active_enrich <- runAme(undiff_active_seq$active_seq,
                               control = "shuffle",
                               database = jaspar_path,
                               outdir = file.path(out_dir, "meme_active/undiff"))

diff_active_seq <- get_active_seq(lm_active_diff, lm_notact_diff, hg.genome)
diff_active_enrich <- runAme(diff_active_seq$active_seq,
                             control = "shuffle",
                             database = jaspar_path,
                             outdir = file.path(out_dir, "meme_active/diff"))

aicar_active_seq <- get_active_seq(lm_active_aicar, lm_nonact_aicar, hg.genome)
aicar_active_enrich <- runAme(aicar_active_seq$active_seq,
                              control = "shuffle",
                              database = jaspar_path,
                              outdir = file.path(out_dir, "meme_active/aicar"))

palm_active_seq <- get_active_seq(lm_active_palm, lm_nonact_palm, hg.genome)
palm_active_enrich <- runAme(palm_active_seq$active_seq,
                             control = "shuffle",
                             database = jaspar_path,
                             outdir = file.path(out_dir, "meme_active/palm"))

# # Allelic sets
# get_allelic_seq <- function(df, df_null, allelic_info, genome) {
#   
#   allelic_coord <- merge(df, allelic_info, by = "refname")
#   
#   allelic_coord <- allelic_coord %>%
#     distinct(rsid, .keep_all = TRUE) %>%
#     dplyr::filter(rsid != "NA") %>%
#     mutate(start = as.numeric(pos) - 100,
#            stop = as.numeric(pos) + 99) %>%
#     distinct(chr, start, stop, .keep_all = TRUE)
#   
#   allelic_coord_gr <- makeGRangesFromDataFrame(allelic_coord)
#   
#   nonall_coord <- merge(df_null, allelic_info, by = "refname")
#   
#   nonall_coord <- nonall_coord %>%
#     distinct(rsid, .keep_all = TRUE) %>%
#     dplyr::filter(rsid != "NA") %>%
#     mutate(start = as.numeric(pos) - 100,
#            stop = as.numeric(pos) + 99) %>%
#     distinct(chr, start, stop, .keep_all = TRUE)
#   
#   nonall_coord_gr <- makeGRangesFromDataFrame(nonall_coord)
#   
#   allelic_seq <- allelic_coord_gr %>%
#     get_sequence(genome)
#   
#   nonallelic_seq = nonall_coord_gr %>%
#     get_sequence(genome)
#   
#   return(list(allelic_seq = allelic_seq,
#               nonallelic_seq = nonallelic_seq))
#   
# }
# 
# undiff_allelic_seq <- get_allelic_seq(lm_allele_undiff, lm_nonall_undiff, allelic_info, hg.genome)
# undiff_allelic_enrich <- runAme(undiff_allelic_seq$allelic_seq,
#                                control = "shuffle",
#                                database = jaspar_path,
#                                outdir = file.path(out_dir, "meme_allelic/undiff"))
# 
# diff_allelic_seq <- get_allelic_seq(lm_allele_diff, lm_nonall_diff, allelic_info, hg.genome)
# diff_allelic_enrich <- runAme(diff_allelic_seq$allelic_seq,
#                              control = "shuffle",
#                              database = jaspar_path,
#                              outdir = file.path(out_dir, "meme_allelic/diff"))
# 
# aicar_allelic_seq <- get_allelic_seq(lm_allele_aicar, lm_nonall_aicar, allelic_info, hg.genome)
# aicar_allelic_enrich <- runAme(aicar_allelic_seq$allelic_seq,
#                               control = "shuffle",
#                               database = jaspar_path,
#                               outdir = file.path(out_dir, "meme_allelic/aicar"))
# 
# palm_allelic_seq <- get_allelic_seq(lm_allele_palm, lm_nonall_palm, allelic_info, hg.genome)
# palm_allelic_enrich <- runAme(palm_allelic_seq$allelic_seq,
#                              control = "shuffle",
#                              database = jaspar_path,
#                              outdir = file.path(out_dir, "meme_allelic/palm"))

###########################
### Plot enrichment res ###
###########################

filter_fp <- function(df) {
  filtered <- df %>%
    filter(!(tp < 2 *fp))
  return(filtered)
}

undiff_active_enrich_jaspar_filt <- filter_fp(undiff_active_enrich_jaspar)
diff_active_enrich_jaspar_filt <- filter_fp(diff_active_enrich_jaspar)
aicar_active_enrich_jaspar_filt <- filter_fp(aicar_active_enrich_jaspar)
palm_active_enrich_jaspar_filt <- filter_fp(palm_active_enrich_jaspar)

colnames(basal_diff_filtered)[8] <- "motif_id"
colnames(aicar_filtered)[8] <- "motif_id"
colnames(palm_filtered)[8] <- "motif_id"

all_active_motifs <- data.frame(c(undiff_active_enrich_jaspar_filt$motif_id,
                                  diff_active_enrich_jaspar_filt$motif_id,
                                  aicar_active_enrich_jaspar_filt$motif_id,
                                  palm_active_enrich_jaspar_filt$motif_id))

undiff_active_enrich_jaspar_filt_summary <- undiff_active_enrich_jaspar_filt %>%
  mutate(log_pval = -log10(pvalue)) %>%
  arrange(desc(log_pval)) %>%
  mutate(label = ifelse(grepl("ZNF", motif_id) | grepl("SP9", motif_id) | grepl("SP3", motif_id),
                        motif_id, paste0("**", motif_id, "**"))) %>%
  mutate(condition = "Undifferentiated")
undiff_active_enrich_jaspar_filt_summary$label <- as.factor(undiff_active_enrich_jaspar_filt_summary$label)
undiff_active_enrich_jaspar_filt_summary$label <- fct_reorder(undiff_active_enrich_jaspar_filt_summary$label,
                                                              undiff_active_enrich_jaspar_filt_summary$pvalue)


diff_active_enrich_jaspar_filt_summary <- diff_active_enrich_jaspar_filt %>%
  mutate(log_pval = -log10(pvalue)) %>%
  arrange(desc(log_pval)) %>%
  mutate(label = ifelse(grepl("ZNF", motif_id) | grepl("SP9", motif_id) | grepl("SP3", motif_id),
                        motif_id, paste0("**", motif_id, "**"))) %>%
  mutate(condition = "Differentiated")
diff_active_enrich_jaspar_filt_summary$label <- as.factor(diff_active_enrich_jaspar_filt_summary$label)
diff_active_enrich_jaspar_filt_summary$label <- fct_reorder(diff_active_enrich_jaspar_filt_summary$label,
                                                            diff_active_enrich_jaspar_filt_summary$pvalue)

aicar_active_enrich_jaspar_filt_summary <- aicar_active_enrich_jaspar_filt %>%
  mutate(log_pval = -log10(pvalue)) %>%
  arrange(desc(log_pval)) %>%
  mutate(label = ifelse(grepl("ZNF", motif_id),
                        motif_id, paste0("**", motif_id, "**"))) %>%
  mutate(condition = "AICAR")
aicar_active_enrich_jaspar_filt_summary$label <- as.factor(aicar_active_enrich_jaspar_filt_summary$label)
aicar_active_enrich_jaspar_filt_summary$label <- fct_reorder(aicar_active_enrich_jaspar_filt_summary$label,
                                                             aicar_active_enrich_jaspar_filt_summary$pvalue)

palm_active_enrich_jaspar_filt_summary <- palm_active_enrich_jaspar_filt %>%
  mutate(log_pval = -log10(pvalue)) %>%
  arrange(desc(log_pval)) %>%
  mutate(label = ifelse(grepl("ZNF", motif_id),
                        motif_id, paste0("**", motif_id, "**")))%>%
  mutate(condition = "Palmitate")
palm_active_enrich_jaspar_filt_summary$label <- as.factor(palm_active_enrich_jaspar_filt_summary$label)
palm_active_enrich_jaspar_filt_summary$label <- fct_reorder(palm_active_enrich_jaspar_filt_summary$label,
                                                            palm_active_enrich_jaspar_filt_summary$pvalue)

# create facet plots
full_active_jaspar_results <- rbind(undiff_active_enrich_jaspar_filt_summary,
                                    diff_active_enrich_jaspar_filt_summary,
                                    aicar_active_enrich_jaspar_filt_summary,
                                    palm_active_enrich_jaspar_filt_summary)

full_active_jaspar_results$condition <- factor(full_active_jaspar_results$condition,
levels = c("Undifferentiated", "Differentiated", "AICAR", "Palmitate"))
facet_colors <- c("Undifferentiated" = "#ff595e", "Differentiated" = "#ffca3a",
                  "AICAR" = "#8ac926", "Palmitate" = "#1982c4")

strip_colors <- strip_themed(background_x = elem_list_rect(fill = c("#ff595e","#ffca3a","#8ac926","#1982c4")))
condition_labels <- c("Undifferentiated" = "Undifferentiated", "Differentiated" = "Differentiated",
                      "AICAR" = "AICAR", "Palmitate" = "Palmitate")

active_enrich_jaspar_filt_plot <- ggplot(full_active_jaspar_results, aes(x = label, y = log_pval)) +
  geom_point(aes(color = tp, fill = tp), size = 4, shape = 21) +
  scale_x_discrete(labels = full_active_jaspar_results$label) +
  geom_point(data = subset(full_active_jaspar_results, adj.pvalue < 0.1),
             aes(x = label, y = log_pval, fill = tp, color = tp),
             size = 4, shape = 21, color = "black", stroke = 1) +
  scale_color_gradientn(colors = c("purple4", "orangered", "yellow1")) +
  scale_fill_gradientn(colors = c("purple4", "orangered", "yellow1")) +
  facet_wrap2(. ~ condition, labeller = labeller(condition = condition_labels),
              strip = strip_colors, ncol = 2, scales = "free_x") +
  labs(x = "Transcription Factor Motif", y = "-log10(p-value)",
       color = "Number of Variants", fill = "Number of Variants") +
  theme_bw(base_family = "Helvetica", base_size = 16) +
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, size = 12),
        strip.background = element_rect(color = NA),
        panel.border = element_blank()) +
  coord_cartesian(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + expand_limits(y = c(-0.1, 23))

active_enrich_jaspar_filt_plot

ggsave(active_enrich_jaspar_filt_plot,
       filename = file.path(fig_dir, "motif_dotplot.png"),
       units = "in", dpi = 600, width = 9, height = 5.4, device = ragg::agg_png())