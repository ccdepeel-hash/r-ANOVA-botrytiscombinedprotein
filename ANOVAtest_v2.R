# 02/2026
# From Anna // Defining paths for argument on farm
args <- commandArgs(trailingOnly = TRUE)

#assign input files to specific variables // adjusting line 6-7 ;) and 21-23
at_combined_file <- args[1]		# Counts file (.csv)
output_dir <- args[2]     # Output directory

# Ensure output directory ends with a slash (for safe concatenation)
if (!grepl("/$", output_dir)) {
	output_dir <- paste0(output_dir, "/")
}

# bringing packages into environment so we can call them
library(tidyr)
library(dplyr)
library(glmmTMB)
library(car)
library(emmeans)

#load count file for arabidopsis
message("Loading BotrytisCombinedProtein file: ", at_combined_file)

#opening and now referring to Dan's proteomic dataset as df; spreadsheet should be in environment (right side) as 'df'.
#at_combined_file vs hardcoding file path; we want to use the argument for the file path instead of hardcoding it in. we can run the script on farm with different files without changing the code
df <- read.csv(file = at_combined_file)

#selecting all data when organism col. is equal to Botroytinia... 
df <- df %>% filter(Organism == "Botryotinia fuckeliana (strain B05.10)")

#From dplyr package; select(!c..) = selecting all EXCEPT ("!c") listed columns
df <- df %>% select(!c(Protein, Protein.Length, Organism, Protein.Existence, Description,X,X.1,X.2,X.3,X.4,Indistinguishable.Proteins))
df <- df %>% pivot_longer(cols = !c(Protein.ID, Gene, Entry.Name), names_to = 'SampleID', values_to= 'Protein_Abundance')
df <- df %>% filter(startsWith(x = SampleID, prefix = 'Col0_B') | startsWith(x = SampleID, prefix = 'tgg12_B'))

#need for anova results
df <- df %>% mutate(Genotype = if_else(condition = startsWith(x = SampleID, prefix = 'Col0'), true = "Col0", false = "tgg12"))

#get a list of proteins we want to get rid of (>=80% 0s for that protein)
provector <- df %>%
  group_by(Protein.ID) %>%
  summarise(percentage = mean(Protein_Abundance == 0), .groups = "drop") %>%
  filter(percentage >= 0.80) %>%
  pull(Protein.ID)
# summarizing to get proportion of samples where genes come out as zero; generated provector

  #filter to remove the low abundance proteins
filterpv <- df %>%
  filter(!Protein.ID %in% provector) #%in% give columns in dataframe and check against some vector, only keep values in this column that show up in established vector

median <- filterpv %>% group_by(SampleID) %>% summarise(Median = median(Protein_Abundance)) #what do we actually want to group our samples by; pipe to summarise. obj creation. need to remove 0s to not skew median 
dfw <- filterpv %>%
  select(!c(Gene, Entry.Name)) %>%
  pivot_wider(names_from = Protein.ID, values_from = Protein_Abundance)
formula <- as.formula("A6RJ45 ~ Genotype")

model <- glmmTMB(data = dfw, formula, family = nbinom2)
# WHERE WE STOPPED
# For-Loop stuff..Get all protein column names.. everything except SampleID and Genotype?
protein_columns <- colnames(dfw)[!colnames(dfw) %in% c("SampleID", "Genotype")]

#creating an empty dataframe to put the looped data into
anova_results <- data.frame(
)

#subset for testing the loop; for loop will run as many times as we tell it to run
# below; protein_columns is whole column.. we are just going to do first five for now. setting that column name to equal only first 5 rows
protein_columns <- protein_columns[1:10]

for(protein in protein_columns) {
  
  # Printing current protein you're in.. for troubleshooting 
  print(paste("Processing:", protein))
  
  #write formula
  formula <- as.formula(paste(protein, "~ Genotype"))
  
  #build model
  model <- glmmTMB(data = dfw, formula, family = nbinom2)
  
  #pull anova results
  anova_result <- Anova(model)
  
  # adding in protein id column with mutate
  anova_result <- anova_result %>% mutate(Protein.ID = protein,
                                          Convergence.Note = model$fit$message)
  
  #append anova result to the full table (hint: use rbind())
  anova_results <- rbind(anova_results, anova_result)
}

# Apply FDR correction to each variable group // "Benjamini-Hochberg (BH) procedure"
anova_results$p_adj <- p.adjust(anova_results$`Pr(>Chisq)`, method = "BH")

# For each protein after model:
# From Anna // Get estimated marginal means on log scale
emm_log <- emmeans(model, ~ Genotype) %>%
  summary() %>%
  as.data.frame() %>%
  dplyr::rename(emmean_log = emmean) %>%
  mutate(Protein.ID = protein)

# From Anna // Get estimated marginal means on response scale
emm_resp <- emmeans(model, ~ Genotype, type = "response") %>%
  summary() %>%
  as.data.frame() %>%
  dplyr::rename(emmean_response = response) %>%
  mutate(Protein.ID = protein)


# From Anna // Do FDR correction (BH)
# anova_fdr <- lapply(anova_split, function(x) {
  # p_values <- x$`Pr(>Chisq)` #pull out p values from your dataframe
  # x$p_adj <- p.adjust(p_values, method = "BH") #making a new column called p_adj
  # return(x)})
# Recombine into single dataframe
# anova_corrected <- do.call(rbind, anova_fdr)
# Reset row names
# rownames(anova_corrected) <- NULL

# bar chart stuff

# Count how many proteins had convergence issues
convergence_summary <- anova_results %>%
  group_by(Convergence.Note) %>%
  summarise(count = n())

# Bar chart
# library(ggplot2)
# ggplot(convergence_summary, aes(x = Convergence.Note, y = count)) +
  # geom_bar(stat = "identity") +
  # labs(title = "Model Convergence Status",
       # x = "Convergence Message",
       # y = "Number of Proteins") +
   # theme_minimal()