library(data.table)
library(dplyr)
library(tidyr)
library(caret)
library(multiROC)
library(PRROC)
library(ggplot2)
library(patchwork)

#--------load data--------
path <- '/home/chorzow/BI/MIKE_SH/bioinf-cov-classifier/epitope-classifier/'

vdjdb <- fread(paste0(path, 'vdjdb_filtered.txt')) %>% 
  filter(species == 'HomoSapiens')  #TODO: ADD FILTERING PROCEDURES FROM THE OTHER SCRIPT

gliph2_output <- fread(paste0(path, 'gliph2/gliph2_output_raw.csv'))
tcrdist3_centroids_output <- fread(paste0(path, 'tcrdist3/tcrdist3_output_raw.csv'))

#--------DATA PREPROCESSING--------
vdjdb$Clonotype <- seq(1:nrow(vdjdb))
#TODO: ADD FILTERING PROCEDURES FROM THE OTHER SCRIPT

#--------GLIPH2--------


gliph2 <- gliph2_output %>% separate(Sample, c('Clonotype', NA), ':') %>% 
  mutate_at('Clonotype', as.numeric) %>%   # we want only numbers from this column
  group_by(index) %>% rename(Cluster = index) %>% 
  mutate(ulTcRb = as.factor(ulTcRb), Scored_class = names(which.max(table(ulTcRb))))  # choose the most prevalent disease in each cluster and assign
#TODO: is this step necessary?
preds <- merge(gliph2, vdjdb, by = 'Clonotype') %>% 
  select(Clonotype, Scored_class, antigen.species, antigen.epitope, Cluster ) %>% 
  rename(Class = antigen.species, Epitope = antigen.epitope) %>% 
  mutate_at(c('Scored_class', 'Class'), as.factor)

metrics <- confusionMatrix(preds$Scored_class, preds$Class)

# assign epitope weights as a share in cluster

weights <- preds %>% group_by(Cluster, Epitope) %>% count(Epitope) %>% 
  group_by(Cluster) %>% 
  summarize(Cluster = Cluster, Epitope = Epitope, Epitope_weight = n / sum(n),
            predicted_epitope = Epitope[which.max(Epitope_weight)])

preds <- merge(weights, preds, by = c('Cluster', 'Epitope'))
preds$correct <- ifelse(preds$Epitope == preds$predicted_epitope, TRUE, FALSE)

roc_gliph <- roc.curve(scores.class0 = preds$Epitope_weight, weights.class0 = preds$correct, curve = TRUE)
pr_gliph <- pr.curve(scores.class0 = preds$Epitope_weight, weights.class0 = preds$correct, curve = TRUE)

roc_gliph_plot <- ggplot(data.frame(roc_gliph$curve), aes(x = X1, y = X2, color = X3)) + 
  geom_line() + scale_color_gradientn(colours = rainbow(10)) + 
  labs(x='FPR', y='Sensitivity', title = 'ROC Curve', 
       subtitle = paste0('AUC = ', format(roc_gliph$auc, digits=3)),
       color = 'Threshold') + 
  theme_bw() + theme(legend.position = 'none')

pr_gliph_plot <- ggplot(data.frame(pr_gliph$curve), aes(x = X1, y = X2, color = X3)) + 
  geom_line() + scale_color_gradientn(colours = rainbow(10)) + 
  labs(x='Precision', y='Recall', title = 'Precision-Recall Curve',
       subtitle = paste0('AUC = ', format(pr_gliph$auc.integral, digits=3)),
       color = 'Threshold') + xlim(0, 1) + ylim(0, 1) + 
  theme_bw()

gliph2_plots <- roc_gliph_plot + pr_gliph_plot + plot_annotation('GLIPH2')

#TODO: add multi-class PR curves with multiROC

#--------TCRDIST3--------
tcrdist3_centroids_output$neighbors <- gsub('\\[|\\]', '', tcrdist3_centroids_output$neighbors)
s <- tcrdist3_centroids_output %>% separate_rows(neighbors, sep = ',')
summ <- s %>% group_by(epitope) %>% 
  summarize(cluster = cluster_id, Clonotype = as.numeric(neighbors))
str(summ)
res <- merge(summ, vdjdb, by = 'Clonotype') %>% 
  select(cluster, epitope, antigen.epitope) %>% 
  rename(predicted_epitope = epitope, epitope = antigen.epitope) %>% 
  group_by(cluster, epitope) %>% mutate(num = n()) %>% 
  group_by(cluster) %>% 
  summarize(cluster = cluster, epitope = epitope, epitope_weight = num / sum(num),
            predicted_epitope = predicted_epitope) %>% arrange(cluster)
res_ded <- res[!duplicated(res),]

res$correct <- ifelse(res$epitope == res$predicted_epitope, TRUE, FALSE)

roc_tcrdist <- roc.curve(scores.class0 = res$epitope_weight, weights.class0 = res$correct, curve = TRUE)
pr_tcrdist <- pr.curve(scores.class0 = res$epitope_weight, weights.class0 = res$correct, curve = TRUE)

roc_tcrdist_plot <- ggplot(data.frame(roc_tcrdist$curve), aes(x = X1, y = X2, color = X3)) + 
  geom_line() + scale_color_gradientn(colours = rainbow(10)) + 
  labs(x='FPR', y='Sensitivity', title = 'ROC Curve', 
       subtitle = paste0('AUC = ', format(roc_tcrdist$auc, digits=3)),
       color = 'Threshold') + 
  theme_bw() + theme(legend.position = 'none')

pr_tcrdist_plot <- ggplot(data.frame(pr_tcrdist$curve), aes(x = X1, y = X2, color = X3)) + 
  geom_line() + scale_color_gradientn(colours = rainbow(10)) + 
  labs(x='Precision', y='Recall', title = 'Precision-Recall Curve',
       subtitle = paste0('AUC = ', format(pr_tcrdist$auc.integral, digits=3)),
       color = 'Threshold') + xlim(0, 1) + ylim(0, 1) + 
  theme_bw()

tcrdist3_plots <- roc_tcrdist_plot + pr_tcrdist_plot + plot_annotation(title = 'Tcrdist3')
gliph2_plots / tcrdist3_plots
