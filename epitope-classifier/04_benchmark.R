library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(pROC)

#--------LOAD DATA--------
path <- '/home/chorzow/BI/MIKE_SH/bioinf-cov-classifier/epitope-classifier/'
source(paste0(path, '00_custom_functions.R'))
source(paste0(path, '01_data_preprocessing.R'))

gliph2 <- fread(paste0(path, 'gliph2/gliph2_output_raw.csv'), fill = T)
tcrdist3 <- fread(paste0(path, 'tcrdist3/centroids_df.csv'))

#--------GLIPH2--------
source(paste0(path, '02_gliph2.R'))

# resulting plots
roc_gliph2_multiclass_plot / pr_gliph2_multiclass_plot  # multiclass
gliph2_multiclass_cm

binary_roc_gliph2_plot / binary_prc_gliph2_plot  # binary
gliph2_binary_cm

#--------TCRDIST3--------
source(paste0(path, '03_tcrdist3.R'))

#resulting plots
roc_tcrdist3_multiclass_plot / pr_tcrdist3_multiclass_plot  # multiclass
tcrdist3_multiclass_cm

roc_tcrdist3_binary_plot / prc_tcrdist3_binary_plot  # binary
tcrdist3_binary_cm
