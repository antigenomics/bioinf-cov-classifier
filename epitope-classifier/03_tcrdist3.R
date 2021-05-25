#--------load data--------
path <- '/home/chorzow/BI/MIKE_SH/bioinf-cov-classifier/epitope-classifier/'
source(paste0(path, '00_custom_functions.R'))
source(paste0(path, '01_data_preprocessing.R'))
tcrdist3 <- fread(paste0(path, 'tcrdist3/centroids_df.csv'))

tcrdist3$neighbors <- gsub('\\[|\\]', '', tcrdist3$neighbors)
summ <- tcrdist3 %>% separate_rows(neighbors, sep = ',') %>% 
  group_by(epitope) %>% summarize(cluster = cluster_id,
                                  clonotype = as.numeric(neighbors))

#--------MAIN--------
res <- merge(summ, vdjdb, by = 'clonotype') %>% 
  select(clonotype, cluster, epitope, antigen.epitope, antigen.species) %>% 
  rename(predicted_epitope = epitope, epitope = antigen.epitope, Species = antigen.species) %>% 
  group_by(cluster, epitope) %>% mutate(num = n()) %>% 
  group_by(cluster) %>% 
  summarize(cluster = cluster, epitope = epitope, epitope_weight = num / sum(num),
            predicted_epitope = predicted_epitope, clonotype = clonotype,
            Species = Species) %>% arrange(cluster) %>%
  mutate(predicted_epitope = epitope[which.max(epitope_weight)]) %>% 
  mutate(ylq_real = ifelse(str_starts(epitope, 'YLQ'), 'YLQ', 'Not_YLQ'),
         ylq_pred = ifelse(str_starts(predicted_epitope, 'YLQ'), 'YLQ', 'Not_YLQ')) %>% 
  mutate_at(c('ylq_real', 'ylq_pred', 'epitope', 'predicted_epitope'), as.factor) %>% 
  mutate(ylq_correct = ifelse(ylq_real == ylq_pred, TRUE, FALSE))


# binary classification: ROC, PRC
tcrdist3_binary_roc <- roc(res$ylq_real, res$epitope_weight)
auc(tcrdist3_binary_roc)
roc_tcrdist3_binary_plot <- ggroc(tcrdist3_binary_roc) + theme_bw() + ggtitle(paste0('tcrdist3 ROC Curve(AUC = ', round(auc(tcrdist3_binary_roc),2), '): Binary classification'))
pr <- coords(tcrdist3_binary_roc, 'all', ret = c('recall', 'precision', 'threshold'), transpose = F)
prc_tcrdist3_binary_plot <- ggplot(pr, aes(recall, precision)) + geom_path(aes(recall, precision)) + coord_equal() + 
  theme_bw() + ggtitle('tcrdist3 PR Curve: Binary classification')

# Multi-class classification: ROC, PRC
tcrdist3_multiclass_roc <- multiclass.roc(res$epitope, res$epitope_weight)
auc(tcrdist3_multiclass_roc)
rs_tcrdist <- tcrdist3_multiclass_roc[['rocs']][c(3, 5, 6)]  # list of roc curves
pr_curves_tcrdist <- lapply(c(1:3), function (i) 
  coords(rs_tcrdist[[i]], 'all', ret = c('recall', 'precision', 'threshold'), transpose = FALSE))  # list of data for PRC

roc_tcrdist3_multiclass_plot <- ggroc(rs_tcrdist) + theme_bw() + scale_color_manual(name="YLQ vs",
                                                                                    labels=c("GIL","GLC","NLP"),
                                                                                    values=c("red","green","blue")) + 
  ggtitle('tcrdist3 ROC Curves: Multi-class classification')

pr_curves_tcrdist[[1]]$response <- 'GIL'
pr_curves_tcrdist[[2]]$response <- 'GLC'
pr_curves_tcrdist[[3]]$response <- 'NLP'

pr_tcrdist3 <- rbind(pr_curves_tcrdist[[1]], pr_curves_tcrdist[[2]], pr_curves_tcrdist[[3]])
pr_tcrdist3_multiclass_plot <- ggplot(pr_tcrdist3, aes(recall, precision, color = response)) + geom_path(aes(recall, precision)) + coord_equal() + 
  theme_bw() + scale_color_manual(name='YLQ vs',
                                  values=c('red', 'green', 'blue')) + 
  ggtitle('tcrdist3 PR Curves: Multi-class classification')

# confusion matrices
tcrdist3_binary_cm <- conf_matrix(res$ylq_real, 
                                  res$ylq_pred, 
                                  title = 'Tcrdist3: Confusion matrix. Binary classification\n')
tcrdist3_multiclass_cm <- conf_matrix(res$epitope, res$predicted_epitope, 
                                      title = 'Tcrdist3 Confusion matrix. Multi-class classification\n')
