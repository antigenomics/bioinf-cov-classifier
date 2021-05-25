#--------load data--------
path <- '/home/chorzow/BI/MIKE_SH/bioinf-cov-classifier/epitope-classifier/'
source(paste0(path, '00_custom_functions.R'))
source(paste0(path, '01_data_preprocessing.R'))

#--------MAIN--------
gliph2 <- fread(paste0(path, 'gliph2/gliph2_output_raw.csv'), fill = T)
tcrdist3 <- fread(paste0(path, 'tcrdist3/centroids_df.csv'))

gliph2 <- gliph2 %>% separate(Sample, c('clonotype', NA), ':') %>%  # extract clonotypes' IDs
  mutate_at('clonotype', as.numeric) %>%   # convert to numeric
  group_by(index) %>% rename(cluster = index) %>%  # group by cluster
  mutate(ulTcRb = as.factor(ulTcRb)) %>% 
  select(-TcRa) %>% filter(!is.na(cluster)) %>% rename(cdr3 = TcRb) %>% 
  filter(number_unique_cdr3 >= 10)

preds <- merge(gliph2, vdjdb, by = c('clonotype', 'cdr3'), all.y = T) %>% # set all.y to T to complete
  select(clonotype, cdr3, antigen.species, antigen.epitope, cluster) %>%  # choose only necessary columns
  rename(species = antigen.species, epitope = antigen.epitope) %>% 
  mutate_at('species', as.factor)
preds$cluster[which(is.na(preds$cluster))] <- 0  # set cluster of unclassified samples to 0

# assign epitope weights as a share in cluster

weights <- preds %>% group_by(cluster, epitope) %>% count(epitope) %>% 
  group_by(cluster) %>% 
  summarize(cluster = cluster, epitope = epitope, epitope_weight = n / sum(n),
            predicted_epitope = epitope[which.max(epitope_weight)])

# attach to the initial dataframe
preds <- merge(weights, preds, by = c('cluster', 'epitope'))
preds$correct_epitope <- ifelse(preds$epitope == preds$predicted_epitope, TRUE, FALSE)

# IMPORTANT OPTION: remove clonotypes that were assigned to multiple clusters
# line with "filter(n == 1) %>%" can be either commented or uncommented depending on whether you want to remove clonotypes that were assigned to multiple clusters (yes, GLIPH2 clusters TCRs this way)

preds <-merge(preds, count(preds, clonotype)) %>% filter(n == 1)

preds$correct_epitope[which(preds$cluster == 0)] <- FALSE
preds$predicted_epitope[which(preds$cluster == 0)] <- 'None'
preds <- preds %>%
  mutate_at(c('epitope', 'predicted_epitope'), as.factor) %>% 
  mutate(ylq_real = ifelse(str_starts(epitope, 'YLQ'), 'YLQ', 'Not_YLQ'),
         ylq_pred = ifelse(str_starts(predicted_epitope, 'YLQ'), 'YLQ', 'Not_YLQ')) %>% 
  mutate_at(c('ylq_real', 'ylq_pred'), as.factor) %>% 
  mutate(ylq_correct = ifelse(ylq_real == ylq_pred, TRUE, FALSE))


# Multi-class classification: ROC, PRC
gliph2_multiclass_roc <- multiclass.roc(preds$epitope, preds$epitope_weight)
auc(gliph2_multiclass_roc)
rs_gliph2 <- gliph2_multiclass_roc[['rocs']][c(3, 5, 6)]  # list of roc curves
pr_curves_gliph2 <- lapply(c(1:3), function (i) 
  coords(rs_gliph2[[i]], 'all', ret = c('recall', 'precision', 'threshold'), transpose = FALSE))  # list of data for PRC

roc_gliph2_multiclass_plot <- ggroc(rs_gliph2) + theme_bw() + scale_color_manual(name="YLQ vs",
                                                                                 labels=c("GIL","GLC","NLP"),
                                                                                 values=c("red","green","blue")) + 
  ggtitle('GLIPH2 ROC Curves: Multi-class classification')

pr_curves_gliph2[[1]]$response <- 'GIL'
pr_curves_gliph2[[2]]$response <- 'GLC'
pr_curves_gliph2[[3]]$response <- 'NLP'

pr_gliph2 <- rbind(pr_curves_gliph2[[1]], pr_curves_gliph2[[2]], pr_curves_gliph2[[3]])
pr_gliph2_multiclass_plot <- ggplot(pr_gliph2, aes(recall, precision, color = response)) + geom_path(aes(recall, precision)) + coord_equal() + 
  theme_bw() + scale_color_manual(name='YLQ vs',
                                  values=c('red', 'green', 'blue')) + 
  ggtitle('GLIPH2 PR Curves: Multi-class classification')

# binary classification: ROC, PRC
gliph2_binary_roc <- roc(preds$ylq_real, preds$epitope_weight)
auc(gliph2_binary_roc)
binary_roc_gliph2_plot <- ggroc(gliph2_binary_roc) + theme_bw() + ggtitle(paste0('GLIPH2 ROC Curve(AUC = ', round(auc(gliph2_binary_roc),2), '): Binary classification'))
pr1 <- coords(gliph2_binary_roc, 'all', ret = c('recall', 'precision', 'threshold'), transpose = F)
binary_prc_gliph2_plot <- ggplot(pr1, aes(recall, precision)) + geom_path(aes(recall, precision)) + coord_equal() + 
  theme_bw() + ggtitle('GLIPH2 PR Curve: Binary classification')

# confusion matrices
gliph2_binary_cm <- conf_matrix(preds$ylq_real, 
                                preds$ylq_pred, 
                                title = 'GLIPH2 Confusion matrix: Binary classification\n')
gliph2_multiclass_cm <- conf_matrix(preds$epitope, 
                                    preds$predicted_epitope, 
                                    title='GLIPH2: Confusion matrix: Multi-class classification\n')
