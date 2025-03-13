# Calculate Heterozygosity from PLINK processed data
# Load the file
GemHumMonHyaRefGemPlink<-read.table("/path/GemHumMonHyaRefGem_geno005_maf001_hwe_header.het", header=T)
#data_Plink処理後
GemHumMonHyaRefGemPlink$He<-c((GemHumMonHyaRefGemPlink$N.NM.-GemHumMonHyaRefGemPlink$E.HOM.)/GemHumMonHyaRefGemPlink$N.NM.)
GemHumMonHyaRefGemPlink$Ho<-c((GemHumMonHyaRefGemPlink$N.NM.-GemHumMonHyaRefGemPlink$O.HOM.)/GemHumMonHyaRefGemPlink$N.NM.)


# Load necessary libraries
library(ggplot2)
library(dplyr)

# Check data
head(GemHumMonHyaRefGemPlink)

# Create boxplots
# Specify the order of species on the x-axis
GemHumMonHyaRefGemPlink$Species <- factor(data$Species, levels = c("hya", "bif", "cyt", "sub", "gem", "hum", "mon", "hyb"))

#Ho
# Plot Ho values
ggplot(GemHumMonHyaRefGemPlink, aes(x = Species, y = Ho)) +
  # Boxplotの中の色を変更
  geom_boxplot(aes(fill = Species), alpha = 0.7, outlier.shape = NA) +
  # データポイントの色を変更
  geom_jitter(aes(color = Species), width = 0.2, size = 1.5, alpha = 0.6) +
  labs(
    title = "Species-wise Boxplot for Ho",
    x = "Species",
    y = "Ho",
    fill = "Species",
    color = "Species"
  ) +
  scale_fill_manual(values = c(
    "hya" = "white",  
    "bif" = "white",  
    "cyt" = "white",  
    "sub" = "white",  
    "gem" = "grey", 
    "hum" = "grey",  
    "mon" = "grey",  
    "hyb" = "grey"   
  )) +
  scale_color_manual(values = c(
    "hya" = "black",  
    "bif" = "black",  
    "cyt" = "black",  
    "sub" = "black",  
    "gem" = "black",  
    "hum" = "black",  
    "mon" = "black",  
    "hyb" = "black"   
  )) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  
    axis.line = element_line(size = 0.5, color = "black"), 
    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.text.y = element_text(size = 10), 
    legend.position = "right" 
  )

#He
ggplot(GemHumMonHyaRefGemPlink, aes(x = Species, y = He)) +
  geom_boxplot(aes(fill = Species), alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = Species), width = 0.2, size = 1.5, alpha = 0.6) +
  labs(
    title = "Species-wise Boxplot for He",
    x = "Species",
    y = "He",
    fill = "Species",
    color = "Species"
  ) +
  scale_fill_manual(values = c(
    "hya" = "white",  
    "bif" = "white",  
    "cyt" = "white",  
    "sub" = "white",  
    "gem" = "grey", 
    "hum" = "grey",  
    "mon" = "grey",  
    "hyb" = "grey"   
  )) +
  scale_color_manual(values = c(
    "hya" = "black",  
    "bif" = "black",  
    "cyt" = "black",  
    "sub" = "black",  
    "gem" = "black",  
    "hum" = "black",  
    "mon" = "black",  
    "hyb" = "black"   
  )) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  
    axis.line = element_line(size = 0.5, color = "black"), 
    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.text.y = element_text(size = 10), 
    legend.position = "right" 
  )

#F
# F box plot
ggplot(GemHumMonHyaRefGemPlink, aes(x = Species, y = F)) +
  geom_boxplot(aes(fill = Species), alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = Species), width = 0.2, size = 1.5, alpha = 0.6) +
  labs(
    title = "Species-wise Boxplot for F",
    x = "Species",
    y = "Inbreefing coefficient (F)",
    fill = "Species",
    color = "Species"
  ) +
  scale_fill_manual(values = c(
    "hya" = "white",  
    "bif" = "white",  
    "cyt" = "white",  
    "sub" = "white",  
    "gem" = "grey", 
    "hum" = "grey",  
    "mon" = "grey",  
    "hyb" = "grey"   
  )) +
  scale_color_manual(values = c(
    "hya" = "black",  
    "bif" = "black",  
    "cyt" = "black",  
    "sub" = "black",  
    "gem" = "black",  
    "hum" = "black",  
    "mon" = "black",  
    "hyb" = "black"   
  )) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  
    axis.line = element_line(size = 0.5, color = "black"), 
    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.text.y = element_text(size = 10), 
    legend.position = "right" 
  )

# Perform random sampling with 4 samples per species and multiple comparisons
# Install necessary packages if not installed
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("multcomp")) install.packages("multcomp")

# load libraries
library(tidyverse)
library(multcomp)
library(emmeans)

# Convert Species column to factor and numeric columns to numeric type
GemHumMonHyaRefGemPlink$Species <- as.factor(GemHumMonHyaRefGemPlink$Species)
GemHumMonHyaRefGemPlink$He <- as.numeric(GemHumMonHyaRefGemPlink$He)
GemHumMonHyaRefGemPlink$Ho <- as.numeric(GemHumMonHyaRefGemPlink$Ho)

#He
# Fit GLM model for He
glm_model <- glm(He ~ Species, data = GemHumMonHyaRefGemPlink)
emmeans_result <- emmeans(glm_model, pairwise ~ Species)
contrast_summary <- as.data.frame(summary(emmeans_result$contrasts))
contrast_summary$p.adj <- p.adjust(contrast_summary$p.value, method = "fdr")
print(contrast_summary)

#Ho
glm_model <- glm(Ho ~ Species, data =GemHumMonHyaRefGemPlink )
emmeans_result <- emmeans(glm_model, pairwise ~ Species)
contrast_summary <- as.data.frame(summary(emmeans_result$contrasts))
contrast_summary$p.adj <- p.adjust(contrast_summary$p.value, method = "fdr")
print(contrast_summary)

#F
glm_model <- glm(F ~ Species, data =GemHumMonHyaRefGemPlink )
emmeans_result <- emmeans(glm_model, pairwise ~ Species)
contrast_summary <- as.data.frame(summary(emmeans_result$contrasts))
contrast_summary$p.adj <- p.adjust(contrast_summary$p.value, method = "fdr")
print(contrast_summary)

##Perform 100 bootstrap iterations for He
if (!require("emmeans")) install.packages("emmeans")
library(emmeans)

results <- replicate(100, {
  tryCatch({
    selected_samples <- GemHumMonHyaRefGemPlink %>%
      group_by(Species) %>%
      sample_n(4)
    
    glm_model <- glm(He ~ Species, data = selected_samples)
    
    tukey_result <- emmeans(glm_model, pairwise ~ Species, adjust = "fdr")
    
    contrast_summary <- summary(tukey_result$contrasts)
    pvalues <- contrast_summary$p.value
    
    return(pvalues)
  }, error = function(e) {
    return(rep(NA, 28))  # Numbers of comparisons
  })
}, simplify = FALSE)

# Convert results to a dataframe
pvalues_df <- do.call(rbind, results) %>%
  as.data.frame()

contrast_summary <- summary(tukey_result$contrasts)
column_names <- contrast_summary$contrast  # 比較名が保存されている列

colnames(pvalues_df) <- column_names

print(head(pvalues_df))

mean_pvalues <- pvalues_df %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))
# 95% CI list
ci_results <- list()

# median and 95CI calculation
for (col_name in colnames(pvalues_df)) {
  ci_values <- quantile(pvalues_df[[col_name]], probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    ci_results[[col_name]] <- c(
    Lower_95CI = ci_values[1],
    Median = ci_values[2],
    Upper_95CI = ci_values[3]
  )
}


ci_df <- as.data.frame(do.call(rbind, ci_results))
ci_df <- cbind(Comparison = rownames(ci_df), ci_df) 
rownames(ci_df) <- NULL 

print(ci_df)
print(mean_pvalues)

#Plot P values
custom_order <- c(
  "hya - cyt", "bif - cyt", "cyt - sub", "hya - bif", "hya - sub", "bif - sub",
  "cyt - gem", "cyt - hum", "cyt - mon", 
  "hya - gem", "hya - hum", "hya - mon", 
  "bif - gem", "bif - hum", "bif - mon", 
  "sub - gem", "sub - hum", "sub - mon","gem - hum","gem - mon","hum - mon",
  "cyt - hyb", "hya - hyb", "bif - hyb", "sub - hyb",
  "gem - hyb", "hum - hyb", "mon - hyb"
)

pvalues_df %>%
  pivot_longer(cols = everything(), names_to = "Comparison", values_to = "Pvalue") %>%
  mutate(Comparison = factor(Comparison, levels = custom_order)) %>% 
  ggplot(aes(x = Comparison, y = Pvalue)) +
  geom_boxplot() +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  ) +
  labs(
    title = "P-value Distribution Across 100 Random Samples",
    x = "Comparison",
    y = "P-value"
  )


# Ho bootstrap
results <- replicate(100, {
  tryCatch({
        selected_samples <- GemHumMonHyaRefGemPlink %>%
      group_by(Species) %>%
      sample_n(4)
    
    glm_model <- glm(Ho ~ Species, data = selected_samples)
    
    tukey_result <- emmeans(glm_model, pairwise ~ Species, adjust = "fdr")
    
    contrast_summary <- summary(tukey_result$contrasts)
    pvalues <- contrast_summary$p.value
    
    return(pvalues)
  }, error = function(e) {
    return(rep(NA, 28)) 
  })
}, simplify = FALSE)

pvalues_df <- do.call(rbind, results) %>%
  as.data.frame()

contrast_summary <- summary(tukey_result$contrasts)
column_names <- contrast_summary$contrast  
colnames(pvalues_df) <- column_names

print(head(pvalues_df))
mean_pvalues <- pvalues_df %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))

ci_results <- list()
for (col_name in colnames(pvalues_df)) {
  ci_values <- quantile(pvalues_df[[col_name]], probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  
  ci_results[[col_name]] <- c(
    Lower_95CI = ci_values[1],
    Median = ci_values[2],
    Upper_95CI = ci_values[3]
  )
}

ci_df <- as.data.frame(do.call(rbind, ci_results))
ci_df <- cbind(Comparison = rownames(ci_df), ci_df) 
rownames(ci_df) <- NULL 

print(ci_df)
print(mean_pvalues)

#Plot bootsrap P of Ho
pvalues_df %>%
  pivot_longer(cols = everything(), names_to = "Comparison", values_to = "Pvalue") %>%
  mutate(Comparison = factor(Comparison, levels = custom_order)) %>%  
  ggplot(aes(x = Comparison, y = Pvalue)) +
  geom_boxplot() +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  ) +
  labs(
    title = "P-value Distribution Across 100 Random Samples",
    x = "Comparison",
    y = "P-value"
  )

#F bootstrap
results <- replicate(100, {
  tryCatch({
       selected_samples <- GemHumMonHyaRefGemPlink %>%
      group_by(Species) %>%
      sample_n(4)
    
    glm_model <- glm(F ~ Species, data = selected_samples)
    
    tukey_result <- emmeans(glm_model, pairwise ~ Species, adjust = "fdr")    
    contrast_summary <- summary(tukey_result$contrasts)
    pvalues <- contrast_summary$p.value
    
    return(pvalues)
  }, error = function(e) {
    return(rep(NA, 28)) 
  })
}, simplify = FALSE)


pvalues_df <- do.call(rbind, results) %>%
  as.data.frame()

contrast_summary <- summary(tukey_result$contrasts)
column_names <- contrast_summary$contrast
colnames(pvalues_df) <- column_names

print(head(pvalues_df))

mean_pvalues <- pvalues_df %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))

print(mean_pvalues)

ci_results <- list()

for (col_name in colnames(pvalues_df)) {
  ci_values <- quantile(pvalues_df[[col_name]], probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  
  ci_results[[col_name]] <- c(
    Lower_95CI = ci_values[1],
    Median = ci_values[2],
    Upper_95CI = ci_values[3]
  )
}

ci_df <- as.data.frame(do.call(rbind, ci_results))
ci_df <- cbind(Comparison = rownames(ci_df), ci_df)  # 比較名を追加
rownames(ci_df) <- NULL  # 行名を削除

print(ci_df)

#Plot p values of bootsrtap analyses of F
pvalues_df %>%
  pivot_longer(cols = everything(), names_to = "Comparison", values_to = "Pvalue") %>%
  mutate(Comparison = factor(Comparison, levels = custom_order)) %>%  # x軸の順序を指定
  ggplot(aes(x = Comparison, y = Pvalue)) +
  geom_boxplot() +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),  # x軸ラベルを45度回転
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  ) +
  labs(
    title = "P-value Distribution Across 100 Random Samples",
    x = "Comparison",
    y = "P-value"
  )

