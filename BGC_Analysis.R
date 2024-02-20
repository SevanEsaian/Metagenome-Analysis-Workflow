#Plotting Genome Size, Abundance, and Diversity --------
library(ggplot2)
BGC.Ab.Div <- read.csv("BGC_Abundance and Diversity.csv")
View(BGC.Ab.Div)
ggplot(BGC.Ab.Div, aes(x = Taxa, y = Total_Length, fill = Taxa)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = c("dodgerblue2", "forestgreen", "yellow3", "orange",
                                "brown", "pink3")) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + ylab("Genome Length (bp)") + xlab("Bacterial Taxa") +
  annotate("text", x=1, y=8500000, label="A") +
  annotate("text", x=2, y=8500000, label="B") +
  annotate("text", x=3, y=8500000, label="C") +
  annotate("text", x=4, y=8500000, label="B") +
  annotate("text", x=5, y=8500000, label="B") +
  annotate("text", x=6, y=8500000, label="D") 
genome.model <- aov(Total_Length~Taxa, data=BGC.Ab.Div)
summary(genome.model)
TukeyHSD(genome.model, conf.level = 95)

ggplot(BGC.Ab.Div, aes(x = Taxa, y = Abundance, fill = Taxa)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = c("dodgerblue2", "forestgreen", "yellow3", "orange",
                               "brown", "pink3")) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + ylab("BGC Abundance per MAG") + xlab("Bacterial Taxa") +
  annotate("text", x=1, y=14, label="A") +
  annotate("text", x=2, y=14, label="B") +
  annotate("text", x=3, y=14, label="C") +
  annotate("text", x=4, y=14, label="A") +
  annotate("text", x=5, y=14, label="B") +
  annotate("text", x=6, y=14, label="D") 
Ab.model <- aov(Abundance~Taxa, data=BGC.Ab.Div)
summary(Ab.model)
TukeyHSD(Ab.model, conf.level = 95)

ggplot(BGC.Ab.Div, aes(x = Taxa, y = Diversity, fill = Taxa)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = c("dodgerblue2", "forestgreen", "yellow3", "orange",
                               "brown", "pink3")) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + ylab("BGC Diversity per MAG") + xlab("Bacterial Taxa") +
  annotate("text", x=1, y=11, label="A") +
  annotate("text", x=2, y=11, label="B") +
  annotate("text", x=3, y=11, label="C") +
  annotate("text", x=4, y=11, label="A") +
  annotate("text", x=5, y=11, label="A") +
  annotate("text", x=6, y=11, label="D") 
Div.model <- aov(Diversity~Taxa, data=BGC.Ab.Div)
summary(Div.model)
TukeyHSD(Div.model, conf.level = 95)

#Abundance and Diversity Regression--------
annotate_text_log <- expression(paste( "R"^2, "=0.91, ", italic("P")<0.001))

ggplot(BGC.Ab.Div, aes(x=Abundance, y=Diversity, color=Taxa))+
  geom_point(size=4, alpha=0.5) +
  scale_color_manual(values = c("dodgerblue2", "forestgreen", "yellow3", "orange",
                               "brown", "pink3")) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size=14),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)) +
  annotate("text", x=1, y=9, label=annotate_text_log) +
  xlab("BGC Abundance per MAG") + ylab("BGC Diversity per MAG")

Ab.Div.reg <- lm(data=BGC.Ab.Div, formula=Abundance~Diversity)
summary(Ab.Div.reg)

#ReDo PCA 10/2/23 -------
RAperBGC <- read.csv("RAperBGC.csv")
View(RAperBGC)

library(tidyr)

pca_data <- RAperBGC %>%
  group_by(Sample, BGC) %>%
  summarize(Relative_Abundance = coalesce(first(Relative_Abundance), 0)) %>%
  pivot_wider(
    id_cols = Sample,
    names_from = BGC,
    values_from = Relative_Abundance,
    values_fill = 0
  )
View(pca_data)

pca_result <- prcomp(pca_data[, -1], scale. = TRUE)

# Extract the PCA scores and add Species_ID as the first column
pca_scores <- as.data.frame(pca_result$x)
pca_scores$Sample <- pca_data$Sample
View(pca_scores)
##write.csv(pca_scores, "PCA_Scores_ReDone_CDS.csv")
# Access the standard deviations (eigenvalues) of each principal component
eigenvalues <- pca_result$sdev^2  # Square the standard deviations

# Calculate the total variance (sum of eigenvalues)
total_variance <- sum(eigenvalues)

# Calculate the proportion of variance explained per axis
variance_explained <- eigenvalues / total_variance

# Calculate the percentage of variance explained
percentage_variance_explained <- (variance_explained * 100)

# Print the percentage of variance explained per axis
print(percentage_variance_explained)


PCA_Scores <- read.csv("PCA_Scores_ReDone_CDS.csv")
View(PCA_Scores)
library(ggplot2)
ggplot(PCA_Scores, aes(x = PC1, y = PC2, color = Taxa, shape = Site)) +
  geom_point(size = 5, alpha = 0.5) +  # Scatter points
  labs(x = "PC1 [28.14 %]", y = "PC2 [20.37 %]") +  # Axis labels
  theme_bw() +
  scale_color_manual(values = c("blue", "forestgreen", "gold2", "darkred", "hotpink2")) +
  theme(
    axis.text = element_text(size = 14),   # Set font size for axis text
    axis.title = element_text(size = 14),  # Set font size for axis titles
    panel.grid.major = element_blank(),    # Remove major grid lines
    panel.grid.minor = element_blank()     # Remove minor grid lines
  ) +
  annotate("text", x = 3.5, y = 4, label = "PERMANOVA (sample): R2 = 0.84, P < 0.001") +
  annotate("text", x = 3.35, y = 3.7, label = "PERMANOVA (BGC): R2 = 0.31, P < 0.001") 
  

# Perform PERMANOVA
View(RAperBGC)
permanova_sample <- adonis(Relative_Abundance ~ Sample, data = RAperBGC)
summary(permanova_sample)
permanova_sample$aov.tab


# Perform PERMANOVA for the 'BGC' factor
permanova_bgc <- adonis(Relative_Abundance ~ BGC, data = RAperBGC)
permanova_bgc$aov.tab

#Cyclobacteriaceae_CDS_Counts-------
Cyclobacteriaceae_cds_counts <- read.csv("Cyclobacteriaceae_cds_counts.csv")
filtered_df <- Cyclobacteriaceae_cds_counts %>%
  select(FullSampleID, BGC)

# Group by FullSampleID and BGC, then summarize to get the sum
summary_df <- filtered_df %>%
  group_by(FullSampleID, BGC) %>%
  summarize(Count_perBGC_perFullSampleID = sum(n(), na.rm = TRUE))

# Spread the summarized data to fill missing BGCs with 0
final_df <- summary_df %>%
  spread(key = BGC, value = Count_perBGC_perFullSampleID, fill = 0)

# Print or save the final dataframe
head(final_df)  # Display the first few rows of the final dataframe
View(final_df)

long_df <- final_df %>%
  pivot_longer(cols = -FullSampleID, names_to = "BGC", values_to = "Count")

# Print or save the long format dataframe
View(long_df)  # Display the first few rows of the long format dataframe

#SRB CDS Counts -----
SRB <- read.csv("Desulfobacterota_SRB_cds_counts.csv")
SRB_filtered_df <- SRB %>%
  select(FullSampleID, BGC)

# Group by FullSampleID and BGC, then summarize to get the sum
SRB_summary_df <- SRB_filtered_df %>%
  group_by(FullSampleID, BGC) %>%
  summarize(Count_perBGC_perFullSampleID = sum(n(), na.rm = TRUE))

# Spread the summarized data to fill missing BGCs with 0
SRB_final_df <- SRB_summary_df %>%
  spread(key = BGC, value = Count_perBGC_perFullSampleID, fill = 0)

# Print or save the final dataframe
  # Display the first few rows of the final dataframe
View(SRB_final_df)

SRB_long_df <- SRB_final_df %>%
  pivot_longer(cols = -FullSampleID, names_to = "BGC", values_to = "Count")

# Print or save the long format dataframe
View(SRB_long_df)  # Display the first few rows of the long format dataframe

#Rhizobiaceae CDS Counts ------
Rhizo <- read.csv("Rhizobiaceae_cds_counts.csv")
Rhizo_filtered_df <- Rhizo %>%
  select(FullSampleID, BGC)

# Group by FullSampleID and BGC, then summarize to get the sum
Rhizo_summary_df <- Rhizo_filtered_df %>%
  group_by(FullSampleID, BGC) %>%
  summarize(Count_perBGC_perFullSampleID = sum(n(), na.rm = TRUE))

# Spread the summarized data to fill missing BGCs with 0
Rhizo_final_df <- Rhizo_summary_df %>%
  spread(key = BGC, value = Count_perBGC_perFullSampleID, fill = 0)

# Print or save the final dataframe
# Display the first few rows of the final dataframe
View(Rhizo_final_df)

Rhizo_long_df <- Rhizo_final_df %>%
  pivot_longer(cols = -FullSampleID, names_to = "BGC", values_to = "Count")

# Print or save the long format dataframe
View(Rhizo_long_df)  # Display the first few rows of the long format dataframe

#Rhodobacteraceae CDS Counts -------
Rhodo <- read.csv("Rhodobacteraceae_cds_counts.csv")
Rhodo_filtered_df <- Rhodo %>%
  select(FullSampleID, BGC)

# Group by FullSampleID and BGC, then summarize to get the sum
Rhodo_summary_df <- Rhodo_filtered_df %>%
  group_by(FullSampleID, BGC) %>%
  summarize(Count_perBGC_perFullSampleID = sum(n(), na.rm = TRUE))

# Spread the summarized data to fill missing BGCs with 0
Rhodo_final_df <- Rhodo_summary_df %>%
  spread(key = BGC, value = Count_perBGC_perFullSampleID, fill = 0)

# Print or save the final dataframe
# Display the first few rows of the final dataframe
View(Rhodo_final_df)

Rhodo_long_df <- Rhodo_final_df %>%
  pivot_longer(cols = -FullSampleID, names_to = "BGC", values_to = "Count")

# Print or save the long format dataframe
View(Rhodo_long_df)  # Display the first few rows of the long format dataframe

#Thiohalocapsa CDS Counts -----
PSB <- read.csv("Thiohalocapsa_PSB_cds_counts.csv")
PSB_filtered_df <- PSB %>%
  select(FullSampleID, BGC)

# Group by FullSampleID and BGC, then summarize to get the sum
PSB_summary_df <- PSB_filtered_df %>%
  group_by(FullSampleID, BGC) %>%
  summarize(Count_perBGC_perFullSampleID = sum(n(), na.rm = TRUE))

# Spread the summarized data to fill missing BGCs with 0
PSB_final_df <- PSB_summary_df %>%
  spread(key = BGC, value = Count_perBGC_perFullSampleID, fill = 0)

# Print or save the final dataframe
# Display the first few rows of the final dataframe
View(PSB_final_df)

PSB_long_df <- PSB_final_df %>%
  pivot_longer(cols = -FullSampleID, names_to = "BGC", values_to = "Count")

# Print or save the long format dataframe
View(PSB_long_df)  # Display the first few rows of the long format dataframe

#Merge them all together ----
merged_df_all_cds_counts <- bind_rows(long_df, SRB_long_df, Rhizo_long_df, Rhodo_long_df, PSB_long_df)
View(merged_df_all_cds_counts)

all_combinations <- expand.grid(
  FullSampleID = unique(merged_df_all_cds_counts$FullSampleID),
  BGC = unique(merged_df_all_cds_counts$BGC)
)
View(all_combinations)

# Merge the original dataframe with the combinations dataframe, filling missing entries with 0
result_df <- left_join(all_combinations, merged_df_all_cds_counts, by = c("FullSampleID", "BGC"))
result_df$Count[is.na(result_df$Count)] <- 0
View(result_df)
All_Merged_Taxa_BGC_Count <- result_df
View(All_Merged_Taxa_BGC_Count)
##write.csv(All_Merged_Taxa_BGC_Count, "Merged BGC Counts for Heatmap.csv")

hm.count <- read.csv("HeatMap_Count_BGC.csv")
View(hm.count)
# Your ggplot code
library(ggplot2)

# Your ggplot code
ggplot(hm.count, aes(x = Taxa, y = BGC, fill = Sample_Count)) +
  geom_tile(color = "black") +
  coord_fixed() +
  geom_text(aes(label = Sample_Count), color = "white", size = 4) +
  scale_fill_gradient(low = "white", high = "darkblue") +
  
  # Customize the theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 11),  # Right-align x-axis labels
    axis.text.y = element_text(size = 14),  # Set font size for y-axis labels
    axis.title = element_text(size = 14)  # Set font size for axis titles
  ) +
  labs(fill="Quantity of MAGs per Taxa")

#Napdos heatmap ------
napdos <- read.csv("NaPDoS Heatmap.csv")
library(ggplot2)
Cyclonapdos<- subset(napdos, Taxa=="Cyclobacteriaceae")
View(Cyclonapdos)
ggplot(napdos, aes(x = Taxa, y = Class_and_BGC, fill = Quantity_of_MAGs)) +
  geom_tile(color = "black") +
  coord_fixed() +
  geom_text(aes(label = Quantity_of_MAGs), color = "white", size = 4) +
  scale_fill_gradient(low = "white", high = "darkblue") +
  theme_bw() +
  # Customize the theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 11),  # Right-align x-axis labels
        axis.text.y = element_text(size = 14),  # Set font size for y-axis labels
        axis.title = element_text(size = 14)  # Set font size for axis titles
  ) +
  labs(fill="Quantity of MAGs per Natural Product") +
  ylab("Natural Product Domain Class")
