Old samples:

- immunity_paper_data - samples from Somuradova et al
-- PBMC - whole blood from COVID+ donors
-- IFNg CD4 - CD4 T-cells that respond to stimulation by coronavirus epitopes
-- IFNg CD8 - CD4 T-cells that respond to stimulation by coronavirus epitopes

New samples (unpublished, should not be shared!):

- convalescent - new data on COVID+ donors (PBMC)
- healthy - healthy donors obtained with the same protocol (PBMC)
- icu - acute COVID+, on intralung ventilation devises (PBMC)
- icu_lung - same as ICU, but T-cells are taken from lung tissue

The idea is to train a classifier that would tell convalescent from healthy, one can use IFNg+ clonotypes (observed in IFNg more than in PBMC) from all data as features for classifier. Other ideas welcome

Can also add data from [Shoukat et al.](https://www.sciencedirect.com/science/article/pii/S2666379121000033)