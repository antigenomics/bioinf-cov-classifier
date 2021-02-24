The goal is to use several TCR specificity classifiers, such as TCRdist and GLIPH to tell A*02 YLQ and A*02 RLQ specific TCRs (those are coronavirus epitopes) from other virus-specific TCRs.

Two algorithms to use for classification are:
1) TCRdist, latest version described in https://www.biorxiv.org/content/10.1101/2020.12.24.424260v1
2) GLIPH, latest version described in https://pubmed.ncbi.nlm.nih.gov/32341563/ 

Here A*02 YLQ means that the data in "mhc.a" column starts with "HLA-A*02" and the "antigen.epitope" column starts with "YLQ" and so on. Note that you should only use TCR beta chain for now: filter data as ``chain == "TRB"``.

As a control set, we will use A*02 GLC, A*02 GIL and A*02 NLV (Influenza, EBV and CMV epitopes respectively).

Please remove datasets with "reference.id" containing "10X" as 10X datasets are too large and we would prefer to draft the classifier performance benchmrk on smaller datasets.

