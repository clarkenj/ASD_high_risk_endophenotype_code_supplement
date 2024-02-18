# Code supplement for the high risk autism phenotype paper
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)
[![DOI](https://img.shields.io/badge/DOI-10.1101%2F2020.06.01.127688%20-informational)](https://doi.org/10.1101/2020.06.01.127688 )

This repository contains the main code used to process and analyse the data presented in the "Reproducible functional connectivity endophenotype confers high risk of ASD diagnosis in a subset of individuals" paper.

## Abstract
Discovery of predictive biomarkers is essential for understanding the neurobiological underpinnings of autism spectrum disorder (ASD) and for improving diagnosis. Biomarkers should confer a high risk of diagnosis at the individual level and be common enough in the general population to be useful. To date, most progress in biomarker detection of ASD has come from the field of genetics, however known genetic risk factors are either common, but associated with a low risk, or associated with a high risk, but extremely uncommon. Functional connectivity (FC) analyses of individuals with ASD have established sensitivity of brain connectivity at the group level. Yet, the translation of these imaging findings into robust markers of individual risk is hampered by the extensive heterogeneity among ASD individuals. Here, we report an FC signature that confers a more than 7-fold increase in individual risk of ASD diagnosis, yet is still identified in an estimated 1 in 200 individuals in the general population (compared to 1 in 90 individuals with ASD). By limiting predictions to only the most confidently identifiable subset of individuals with ASD we were able to increase the individual risk of our prediction by more than 3-fold over that of previously published imaging models. The identified FC risk signature was characterised by underconnectivity of transmodal brain networks and generalised to independent data. Our results demonstrate the ability of a highly targeted prediction model to meaningfully decompose part of the heterogeneity of the autism spectrum. The identified FC signature may help better delineate the multitude of etiological pathways and behavioural symptoms that challenge our understanding of the autism spectrum.

![Figure 2](https://github.com/clarkenj/ASD_high_risk_endophenotype_code_supplement/blob/master/fig2_model.png)

**Figure 2**. HRS is more common than genetic risk markers, confers higher risk than traditional imaging models, and meets the current state-of-the-art in machine learning. a) Monogenic syndromes (green rhombs) and recurrent Copy Number Variants (pink triangles) confer high risk of ASD diagnosis (vertical axis), but are rare (horizontal axis). ASD related single nucleotide polymorphisms (yellow triangles) are very common, but confer negligible risk of ASD. Current imaging based predictive models (pink circles) identify large portions of the general population with low risk of ASD. The high risk ASD signature (orange, black outline) identifies a small portion of the general population with elevated risk of ASD diagnosis, concordant with the estimated performance in the discovery data (orange plus signs), meeting the PPV of 10 machine learning models combined (red circle), using a simple model.

## Scripts
- `Discovery_Conformal_Score.R` and `Validation_Conformal_Score.R` generate conformal scores for ASD and NTC individuals for each network, using bootstrapping, on the training (ABIDE 1) and validation (ABIDE 2) datasets respectively. Due to the large number of files generated by results on the training data, we share results for one of 100 bootstrap iterations as an example, in [Data](./Data). All results for the validation can be found in [20190524_Validation_Data](./20190524_Validation_Data).

- `Discovery_Read_Conformal_Scores.R` and `Validation_Read_Conformal_Scores.R` read the bootstrapped scores and compute combined scores (see paper for details). Results can be found in the respective directories.