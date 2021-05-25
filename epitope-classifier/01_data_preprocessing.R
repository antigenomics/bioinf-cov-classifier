library(dplyr)
library(stringr)


path <- '/home/chorzow/BI/MIKE_SH/bioinf-cov-classifier/epitope-classifier/'

#--------VDJDB FILTERING--------

vdjdb <- fread(paste0(path, 'vdjdb.slim.txt'))

vdjdb <- vdjdb[-grep('10x', vdjdb$reference.id),] %>% 
  filter(gene == 'TRB', species == 'HomoSapiens') %>% filter(
         (str_starts(antigen.epitope, 'YLQ') & str_starts(mhc.a, 'HLA-A\\*02'))
         | (str_starts(antigen.epitope, 'GLC') & str_starts(mhc.a, 'HLA-A\\*02'))
         | (str_starts(antigen.epitope, 'GIL') & str_starts(mhc.a, 'HLA-A\\*02'))
         | (str_starts(antigen.epitope, 'NLV') & str_starts(mhc.a, 'HLA-A\\*02')))
vdjdb$clonotype <- seq(1:nrow(vdjdb))

#--------TCRDIST3--------

tcrdist3_input <- vdjdb %>% 
  mutate(cohort = ifelse(antigen.species == 'SARS-CoV-2', 'CoV', 'non-CoV'),
         count = 1) %>% rename(epitope = antigen.epitope,
                               v_b_gene = v.segm,
                               j_b_gene = j.segm,
                               cdr3_b_aa = cdr3,
                               subject = clonotype) %>% 
  select(subject, epitope, count, v_b_gene, j_b_gene, cdr3_b_aa)
write.table(tcrdist3_input, paste0(path, 'tcrdist3/tcrdist3_input.tsv'), sep = '\t', 
            quote = F, row.names = F, col.names = T)
rm(tcrdist3_input)
# tcrdist3 is run using the attached script tcrdist3_clustering.py. 
# see tcrdist3_clustering.ipynb for details and intermediate results. 

#--------GLIPH2--------
gliph2_input <- vdjdb %>% 
  mutate(cohort = ifelse(antigen.species == 'SARS-CoV-2', 'CoV', 'non-CoV'),
         count = 1) %>% rename(CDR3b = cdr3,
                               TRBV = v.segm,
                               TRBJ = j.segm,
         ) %>% unite('subject:condition', clonotype:cohort, sep=':') %>%
  mutate(count = 1, CDR3a = NA) %>% 
  select(CDR3b, TRBV, TRBJ, CDR3a, 'subject:condition', count) 

write.table(gliph2_input, paste0(path, 'gliph2/gliph2_input.tsv'), sep = '\t',
            quote = F, row.names = F, col.names = T)
rm(gliph2_input)
# run gliph2: via web interface @ http://50.255.35.37:8080

#--------GLIPH--------
# gliph_input <- vdjdb %>% select(gene, cdr3, v.segm, j.segm, antigen.epitope) %>% 
#   rename(CDR3b = cdr3, TRBV = v.segm, TRBJ = j.segm, epitope = antigen.epitope)
# 
# write.table(gliph_input[, c(2:4)], file = paste0(path, 'gliph_input.txt'), 
#             quote = F, row.names = F, sep = '\t')

# run gliph: Software/gliph/bin/gliph-group-discovery.pl --tcr bioinf-cov-classifier/epitope-classifier/TCRB_gene.txt
