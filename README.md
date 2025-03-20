# Mycobiome-Bacteriome Analysis

## Project Overview
This repository contains scripts and data for analyzing the gut mycobiome-bacteriome interface across humans and nonhuman primates with different subsistence strategies. The fungal composition was assessed using ITS2 sequencing, while the bacterial composition was determined using the V4 region of the 16S rRNA gene. Standard Qiime2 pipelines were used for processing after adapter removal and quality filtering.

## Data Processing Workflow
### 1. Quality Control & Adapter Removal
- Raw sequencing reads were quality-checked using `FastQC`.
- Adapters and low-quality sequences were removed using `cutadapt`.

### 2. Sequence Processing Using Qiime2
#### Import Data into Qiime2:
```bash
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest_file.tsv \
  --output-path demux.qza \
  --input-format PairedEndFastqManifestPhred33V2
```
#### Denoising with DADA2:
```bash
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left-f 10 --p-trim-left-r 10 \
  --p-trunc-len-f 250 --p-trunc-len-r 250 \
  --o-table feature-table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza
```
#### Taxonomic Assignment:
For bacterial sequences (16S rRNA):
```bash
qiime feature-classifier classify-sklearn \
  --i-classifier greengenes_classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza
```
For fungal sequences (ITS2):
```bash
qiime feature-classifier classify-sklearn \
  --i-classifier unite_classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza
```
#### Phylogenetic Tree Construction:
```bash
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
```

## Repository Structure
```
mycobiome_bacteriome/
│-- README.md  # Project documentation
│-- scripts/   # Additional scripts for data processing
│   │-- BAC_Taxa-analysis.R
│   │-- Fungi_Taxonomy_plots.R
│   │-- Bubble_ITS.R
│   │-- Distance_from_Centroid.R
│
│── data/
│   ├── 16S_feature-table.txt
│   ├── 16S_taxonomy.tsv
│   ├── 16S_tree.nwk
│   ├── ITS_feature-table.txt
│   ├── ITS_taxonomy.tsv
│   ├── ITS_taxonomy_formatted.tsv
│   ├── ITS_tree.nwk
│   ├── Metadata_Filtered.txt
│   ├── Metadata_Filtered_Only_for_Phyloseq_function.txt
│
│── processed_data/
│   ├── 16S_ASV_filtered_Rel_Abun.txt
│   ├── 16S_Bray_distance-matrix_1kqiime.tsv
│   ├── Bray_Clusters_Group_Counts.xlsx
│   ├── Cumulative_ITS_rel_abun.txt
│   ├── Cumulative_combined_relabun.txt
│   ├── CytoScapeStats-For_AllGroups.txt
│   ├── ITS_ASV_filtered_Rel_Abun.txt
│   ├── ITS_Bray_distance-matrix_1kqiime.tsv
│   ├── ITS_family_proportions.txt
│   ├── Sel_16S_ASV_filtered_Rel_Abun.txt
│   ├── Sel_30_genus_rel_foramtted.txt
│   ├── Sel_ITS_ASV_filtered_Rel_Abun.txt
│   ├── Total_Diversity.txt
│   ├── procrustes.qzv
│
│── scripts/
│   ├── BAC_Taxa-analysis.R
│   ├── BaAka_all_Corr.R
│   ├── Bacteria_Taxonomy_plots.R
│   ├── Bubble_16S.R
│   ├── Bubble_ITS.R
│   ├── Bubble_ITS_reviewer.R
│   ├── Cohesion_From_matrix.R
│   ├── Distance_from_Centroid.R
│   ├── Fungi_Taxonomy_plots.R
│   ├── ITS_Taxa-analysis.R
│   ├── Network_Stats.R
│   ├── Network_attributes.R
│   ├── PCoA_ITS_after_Removing_Gorilla_repeats.R
│   ├── Procrustes.R
│   ├── Revised_Figures_removeCpativeH.R
│   ├── Taxa-analysis.R
│   ├── pick_details.pl
│
│── figures/
│   ├── Figure1-2_Alpha-Beta-diversity.R
│   ├── Figure3A_and5.R
│   ├── Figure3B.R
│   ├── Figure4_new.R
```

The following files were used for downstream statistical analysis
1. Bacterial and fungal ASV abundances
2. Tree file (rooted)
3. Taxonomy file (using UNITE database for Fungal and GreenGenes for bacterial)
4. Metadata file

## Downstream Statistical Analysis
The following R scripts were used for generating figures and statistical analyses:
- **Alpha and Beta Diversity**: `Figure1-2_Alpha-Beta-diversity.R`
- **Bacterial Fungal correlation**: `Figure3A_and5.R`
- **Network analysis and statistics**: `Figure3B.R` and `Figure4_new.R`
- **Taxonomy Plots**:
  - Bacteria: `Bacteria_Taxonomy_plots.R`
  - Fungi: `Fungi_Taxonomy_plots.R`
- **Network Analysis**:
  - `Network_Stats.R`
  - `Network_attributes.R`
For detailed analysis workflows and figure generation, refer to the scripts in the `scripts/` directory.

## References
- Qiime2 Documentation: [https://docs.qiime2.org](https://docs.qiime2.org)
- For Databases: [https://docs.qiime2.org/2024.10/data-resources/](https://docs.qiime2.org/2024.10/data-resources/) 
- GreenGenes Database for Bacterial Taxonomy: [https://ftp.microbio.me/greengenes_release](https://ftp.microbio.me/greengenes_release)
- UNITE Database for Fungal Taxonomy: [https://unite.ut.ee](https://unite.ut.ee)

## Citation
If you use this repository, please cite the original research article.
If you want to access the publication, here is the title and link: [https://www.nature.com/articles/s41522-022-00274-3](https://www.nature.com/articles/s41522-022-00274-3)

## Contact
For any queries, please reach out via compbiosharma@gmail.com OR ashoks773@gmail.com

---
This README provides an overview of the dataset and analysis pipeline. For further details, refer to individual scripts and data files.
