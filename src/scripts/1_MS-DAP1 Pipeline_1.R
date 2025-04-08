library(msdap)
############## Starting msdap and loading raw files##################
dataset = import_dataset_diann(paste0(diann_data))
#dataset = import_dataset_spectronaut(paste0(spectronautdata))
dataset = import_fasta(dataset,files = c(paste0(FASTA)))
#write_template_for_sample_metadata(dataset, paste0(Metadata))
dataset = import_sample_metadata(dataset, filename = paste0(Metadata))
dataset = setup_contrasts(
  dataset,
  contrast_list = list(
    c(paste0(NameCond2), paste0(NameCond1)),
    c(paste0(NameCond1), paste0(NameCond2))
  )
)

                                                          #c("NC", "Htt"), 
                                                          #c("NC", "Wt")


############## Running MS-DAP Pipeline Part 1####################
my_norm = function(x_as_log2, ...) {
  cat("Example custom normalization implementation, my_norm(), returning log2 peptides without any norm")
  return(x_as_log2)
}
# 6) Main function that runs the entire pipeline
dataset = analysis_quickstart(
  dataset,
  filter_min_detect = 1,
  filter_min_quant = 1,
  filter_min_peptide_per_prot = 1, #protein-level filter; at least 1 peptide
  norm_algorithm = "my_norm", # normalization using custom algorithm = no norm here
  dea_algorithm = c("msqrob"), # statistics; apply multiple methods in parallel/independently
  dea_qvalue_threshold = 0.05,                      # threshold for significance of adjusted p-values in figures and output tables
  dea_log2foldchange_threshold = 0,                # threshold for significance of log2 foldchanges. 0 = disable, NA = automatically infer through bootstrapping
  output_qc_report = TRUE,                          # optionally, set to FALSE to skip the QC report (not recommended for first-time use)
  output_abundance_tables = TRUE,                   # optionally, set to FALSE to skip the peptide- and protein-abundance table output files
  output_dir = paste0(FolderName), # output directory, here set to "msdap_results" within your working directory. Alternatively provide a full path, eg; output_dir="C:/path/to/myproject"
  output_within_timestamped_subdirectory = TRUE )

# print a short summary of results at the end
print_dataset_summary(dataset)
print("done")
# print a short summary of results at the end
print_dataset_summary(dataset)

#  ######## ELLA PARAMETERS #########
# dataset = analysis_quickstart(
#   dataset,
#   filter_min_detect = 1,            # each peptide must have a good confidence score in at least N samples per group
#   filter_min_quant = 1,             # similarly, the number of reps where the peptide must have a quantitative value
#   #filter_fraction_detect = 0.1,    # each peptide must have a good confidence score in at least 75% of samples per group
#   #filter_fraction_quant = 0.1,     # analogous for quantitative values
#   filter_min_peptide_per_prot = 2,  # minimum number of peptides per protein (after applying above filters) required for DEA. Set this to 2 to increase robustness, but note that'll discard approximately 25% of proteins in typical datasets (i.e. that many proteins are only quantified by 1 peptide)
#   filter_by_contrast = TRUE,        # only relevant if dataset has 3+ groups. For DEA at each contrast, filters and normalization are applied on the subset of relevant samples within the contrast for efficiency, see further MS-DAP manuscript. Set to FALSE to disable and use traditional "global filtering" (filters are applied to all sample groups, same data table used in all statistics)
#   norm_algorithm = c(""), # normalization; first vsn, then modebetween on protein-level (applied sequentially so the MS-DAP modebetween algorithm corrects scaling/balance between-sample-groups). "mwmb" is a good alternative for "vsn"
#   dea_algorithm = c("deqms"), #"msempire", "msqrob"), # statistics; apply multiple methods in parallel/independently
#   dea_qvalue_threshold = 0.05,                      # threshold for significance of adjusted p-values in figures and output tables
#   dea_log2foldchange_threshold = 0,                # threshold for significance of log2 foldchanges. 0 = disable, NA = automatically infer through bootstrapping
#   #diffdetect_min_peptides_observed = 2,             # 'differential detection' only for proteins with at least 2 peptides. The differential detection metric is a niche usecase and mostly serves to identify proteins identified near-exclusively in 1 sample group and not the other
#   #diffdetect_min_samples_observed = 3,              # 'differential detection' only for proteins observed in at least 3 samples per group
#   #diffdetect_min_fraction_observed = 0.5,           # 'differential detection' only for proteins observed in 50% of samples per group
#   output_qc_report = TRUE,                          # optionally, set to FALSE to skip the QC report (not recommended for first-time use)
#   output_abundance_tables = TRUE,                   # optionally, set to FALSE to skip the peptide- and protein-abundance table output files
#   output_dir = paste0(FolderName),                     # output directory, here set to "msdap_results" within your working directory. Alternatively provide a full path, eg; output_dir="C:/path/to/myproject",
#   output_within_timestamped_subdirectory = TRUE
# )

