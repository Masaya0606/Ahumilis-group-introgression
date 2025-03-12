library(FactoMineR)
library(factoextra)
library(vegan)

# Load data
eigenvec_file <- "~/NGS_analyses/plink_mac_20220402/GemHumMon2401geno005maf001Hwe/GemHumMon2401refAgemv_geno005_maf001_hwe_permanova.eigenvec"
eigenval_file <- "~/NGS_analyses/plink_mac_20220402/GemHumMon2401geno005maf001Hwe/GemHumMon2401refAgemv_geno005_maf001_hwe_ed.eigenval"

df <- read.table(eigenvec_file, header = TRUE, stringsAsFactors = FALSE)
colnames(df) <- c("FID", "IID", paste0("PC", 1:20))

# Load eigenvalues
eigenvalues <- read.table(eigenval_file, header = FALSE)
eigenvalues <- eigenvalues$V1  # Eigenvalue vector

# Extract group information
groups <- df$FID

# Extract principal component scores
pca_scores <- df[, 3:21]

# Compute cumulative variance
cumulative_var <- cumsum(eigenvalues) / sum(eigenvalues)

# Select the number of principal components that explain at least 70% of the variance
selected_pcs <- which(cumulative_var >= 0.7)[1]
print(selected_pcs)

# Extract selected principal component scores
selected_scores <- pca_scores[, 1:selected_pcs]

# Create a dataframe
df_selected <- as.data.frame(selected_scores)
df_selected$group <- groups

# Perform MANOVA
manova_result <- manova(as.matrix(df_selected[, -ncol(df_selected)]) ~ group, data=df_selected)

# Display MANOVA results
summary(manova_result)

# Perform PERMANOVA
permanova_result <- adonis2(as.matrix(df_selected[, -ncol(df_selected)]) ~ group, data=df_selected, method="euclidean", permutations=999)

# Display PERMANOVA results
print(permanova_result)

# Obtain group pairs
group_levels <- unique(df_selected$group)
pairwise_results <- list()
p_values <- c()  # Vector to store p-values

# Pairwise comparison between groups
for (i in 1:(length(group_levels) - 1)) {
  for (j in (i + 1):length(group_levels)) {
    group_i <- group_levels[i]
    group_j <- group_levels[j]
    
    # Filter data for groups i and j
    df_pair <- df_selected[df_selected$group %in% c(group_i, group_j), ]
    
    # Perform PERMANOVA
    result <- adonis2(as.matrix(df_pair[, -ncol(df_pair)]) ~ group, data=df_pair, method="euclidean", permutations=min(999, choose(nrow(df_pair), 2)))
    
    # Store results
    pairwise_results[[paste(group_i, "vs", group_j)]] <- result
    p_values <- c(p_values, result$aov.tab$`Pr(>F)`[1])  # Store p-values
  }
}

library(vegan)

# Control permutation in PERMANOVA
set.seed(123)  # Ensure reproducibility
ctrl <- how(nperm = 999)  # Control permutations

# Pairwise comparisons between groups
for (i in 1:(length(group_levels) - 1)) {
  for (j in (i + 1):length(group_levels)) {
    group_i <- group_levels[i]
    group_j <- group_levels[j]
    
    # Filter data for groups i and j
    df_pair <- df_selected[df_selected$group %in% c(group_i, group_j), ]
    
    # Perform PERMANOVA (with permutation control)
    result <- adonis2(as.matrix(df_pair[, -ncol(df_pair)]) ~ group, 
                      data=df_pair, 
                      method="euclidean", 
                      permutations=ctrl)
    
    # Store results
    pairwise_results[[paste(group_i, "vs", group_j)]] <- result
    p_values <- c(p_values, result$`Pr(>F)`[1])
  }
}

library(pairwiseAdonis)

# FDR correction
pairwise.adonis2(as.matrix(df_selected[, -ncol(df_selected)]) ~ group, data = df_selected, permutations = 999, method = "euclidean", p.adjust.m = "fdr")
