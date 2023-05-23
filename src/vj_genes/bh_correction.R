
# Create new directory
create_directory <- function(directory_name){
  if (!dir.exists(directory_name)) {
    dir.create(directory_name)
    print(paste0('Create ',directory_name))
  }
}


# Results directory
corrected_dir <- here::here('results','vj_genes','BH_corrected_results')
create_directory(corrected_dir)


# Multiple testing correction
multiple_testing_correction <- function(results){
    results[,'adjusted_p_value'] <- p.adjust(results[,'p_value'], method = 'BH')
  return(results)
}


# Save results
save_results <- function(results, results_dir, file_name){
  print(paste0('Save results in: ',results_dir,'/',file_name))
  write.csv(results,paste0(results_dir,'/',file_name))
}


# Perform multiple testing correction for all WT1 genes
filenames <- c('WT1-37_enriched_TRBJ_gene.tsv','WT1-37_enriched_TRBV_gene.tsv','WT1-126_enriched_TRBJ_gene.tsv','WT1-126_enriched_TRBV_gene.tsv')

for (fn in filenames) {
  # Read data
  results <- read.csv(here::here('results','vj_genes',fn),sep=',',na = 'NA')
  # Perform multiple testing correction
  corrected_results <- multiple_testing_correction(results)
  # Save results
  save_results(corrected_results,corrected_dir, fn)
}


