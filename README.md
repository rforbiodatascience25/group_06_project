# group_06_project

## Project Contributors

**Alberte Englund**

student-id: s215067

github-username: albedamm

**Mathilde Due**

student-id: s215063

github-username: MathildeDue

**Line Winther Gormsen**

student-id: s215111

github-username: LineWintherGormsen

**Sigrid Frandsen**

student-id: s205875

github-username: sigridfrandsen

**Kristine Johansen**

student-id: s215098

github-username: s215098

## Presentation link

Here is a direct link to our presentation:

[`https://raw.githack.com/rforbiodatascience25/group_06_project/main/doc/presentation.html`](https://raw.githack.com/rforbiodatascience25/group_06_project/main/doc/presentation.html)

## Background

The human vaginal microbiome plays an essential role in reproductive health and protection against infections. It is usually dominated by *Lactobacillus* species. However, the vaginal microbiome can differ between people and across health states. Other bacteria, such as Gardnerella, Atopobium, and Prevotella, often appear during dysbiosis or bacterial vaginosis. [1, 2]

Recent advances in metagenomic sequencing enable researchers to reconstruct individual microbial genomes from complex communities. These reconstructed genomes are called metagenome-assembled genomes (MAGs). MAGs provide new insights into the genetic makeup and functions of microorganisms that cannot be easily grown in the lab. Databases such as MGnify collect MAGs from many human and environmental microbiomes, including the vaginal microbiome, and provide useful but often unstandardized metadata.

The dataset used in this project, called `genomes-all_metadata.tsv`, comes from MGnify’s vaginal microbiome genome catalogue. It contains thousands of bacterial MAGs with information on assembly quality (genome length, N50, GC content, completeness, contamination), taxonomic classification, and limited geographical data. Like many real biological datasets, it includes inconsistencies, missing values, and mixed formats [3].

[1] MacSharry, J., Kovács, Z., Xie, Y. *et al.* Endometriosis specific vaginal microbiota links to urine and serum *N*-glycome. *Sci Rep* **14**, 25372 (2024). <https://doi.org/10.1038/s41598-024-76125-2>.

[2] Ravel, J., Gajer, P., Abdo, Z. et al. Vaginal microbiome of reproductive-age women. *Proc Natl Acad Sci USA* 108 (Suppl 1), 4680–4687 (2011). <https://doi.org/10.1073/pnas.1002611107>.

[3] MGnify. Human Vaginal Microbiome Genome Catalogue v1.0. European Bioinformatics Institute (EMBL-EBI). (2023). <https://www.ebi.ac.uk/metagenomics/genome-catalogues/human-vaginal-v1-0/>

### Aim

The aim of this project is to clean, explore, and analyse the genome meta-data-set from the human vaginal microbiome. With this project we aim to uncover patterns in genome quality, taxonomic composition, and ecological characteristics, while demonstrating the principles of reproducible and collaborative data science in R.

In this project we will:

1.  Tidy and standardize the raw meta-data-set by addressing missing values, harmonizing column formats, and parsing the hierarchical taxonomic lineage into separate levels (phylum to species).
2.  Describe the cleaned and augmented dataset, including basic summaries of genome quality, taxonomy, and geography.
3.  Perform exploratory analyses of genome quality and taxonomic composition
4.  Analyse whether endometriosis-associated MAGs are genomically distinct from non-associated MAGs
5.  Investigate geographical patterns in taxonomic composition and endometriosis association across countries using heatmaps, PCA and statistical tests.

Through these analyses, we aim to gain a deeper understanding of the genomic diversity of the vaginal microbiome, while showcasing the entire data science workflow using the tidyverse systems in R.

## Repository Structure

#### data/

-   raw/
    -   raw_data.xlsx
-   01_dat_load.tsv
-   02_dat_clean.tsv
-   03_dat_aug.tsv

#### R/

-   00_all.qmd
-   01_load.qmd
-   02_clean.qmd
-   03_augment.qmd
-   04_describe.qmd
-   05_analysis_1.qmd
-   06_analysis_2.qmd
-   07_analysis_3.qmd
-   08_analysis_4.qmd
-   99_proj_func.R

#### results/

-   \*.html
-   04_key_plot_1.png
-   04_key_plot_2.png
-   04_key_plot_3.png
-   05_key_plot_1.png
-   06_key_plot_1.png
-   07_key_plot_1.png
-   07_key_plot_2.png
-   08_key_plot_1.png
-   08_key_plot_2.png
-   08_key_plot_3.png
-   08_key_plot_4.png
-   08_key_plot_5.png
-   08_key_table_1.png

#### doc/

-   presentation.qmd
-   presentation.html

## Data retrieval

The data for this project is publicly available from the EBI MGnify database.

Our 01_load.qmd script downloads the metadata file automatically and saves it to data/\_raw/.

Direct link to the dataset:

<https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-vaginal/v1.0/genomes-all_metadata.tsv>

## Additional data

During augmentation of the data we are annotating MAGs with endometriosis-associated genera based on the article:

[1] MacSharry, J., Kovács, Z., Xie, Y. *et al.* Endometriosis specific vaginal microbiota links to urine and serum *N*-glycome. *Sci Rep* **14**, 25372 (2024). <https://doi.org/10.1038/s41598-024-76125-2>.
