# Glioma-IDH-state
here , I reproduced this paper https://www.cell.com/cell/fulltext/S0092-8674(15)01692-X   for unsupervised clustering of gliomas (n = 516 LGG) based on methylation levels to classify the six IDH statuses as reported in the paper.
**Glioma Classification Using kNN and Methylation Clustering**

 **1\. Introduction to Gliomas, IDH Status, and Methylation**

Gliomas are aggressive brain tumors categorized by their IDH (Isocitrate Dehydrogenase) mutation status, which influences prognosis and treatment outcomes \[1\]. IDH-mutant gliomas generally have better clinical outcomes. DNA methylation, an epigenetic modification, plays a critical role in gene expression and tumor behavior, making it a useful biomarker for glioma classification \[2\]. Analysis of methylation patterns helps identify tumor subtypes and predict patient prognosis \[3\].

 **2\. Dataset and Preprocessing**

We utilized DNA methylation data from lower-grade gliomas (LGG) available from The Cancer Genome Atlas (TCGA), collected via the Illumina Human Methylation 450 platform. The dataset included metadata on IDH status and glioma clusters.

 Preprocessing Steps:

1\. Removal of missing values.

2\. Exclusion of samples without IDH status.

3\. Application of a mean beta value threshold of 0.3 to filter significant methylation sites.

4\. Merging methylation data with clinical metadata, including IDH status and clustering information.

 **3\. Methodology for KNN Implementation**

To classify the glioma subtypes, we applied the k-Nearest Neighbors (kNN) algorithm, a distance-based classification technique.

 Steps:

1\. Data Preparation: 

   \- The dataset was merged with metadata, and non-clustering columns were removed.

   \- Standardization was applied, transforming the data to have a mean of 0 and a standard deviation of 1\.

2\. Feature Selection: 

   \- Non-relevant columns (e.g., IDH status) were excluded.

3\. Elbow Method: 

   \- This was used to determine the optimal number of clusters by plotting inertia against k values (see Figure 1).
   <img width="602" height="388" alt="elbow" src="https://github.com/user-attachments/assets/53945003-9cab-429f-8e38-9c4b797e0ae8" />
   
Fig. 1: Elbow Method indicates the optimal number of clusters 

4\. Clustering:

   \- kNN was applied, and k=3 was selected for clustering based on the elbow plot results.

5\. Dimensionality Reduction: 

   \- Principal Component Analysis (PCA) was used to visualize the clustering patterns and reduce dimensionality.

 **4\. Results of kNN**

 Clustering Outcome:

\- Three Clusters Identified: 

   \- PCA visualized the clusters based on IDH status (mutant and wild-type), showing that IDH mutations are the primary drivers of variance in the data (see Figure 2). However, some overlap remained, likely due to glioma subtype variability.
   <img width="627" height="508" alt="knn" src="https://github.com/user-attachments/assets/1d572901-f902-409f-b28c-2f28095ee397" />
   
   Fig. 2: PCA shows clusters for IDH mutant 

**\- Cluster Interpretation:** 

   \- The clustering showed a split between IDH-mutant and IDH wild-type groups, aligning with known biological distinctions in gliomas.
<img width="715" height="278" alt="cluster" src="https://github.com/user-attachments/assets/7ca4be94-8b5b-4b09-8849-5d2ede2f347e" />

Fig. 3: Cluster characteristics by grouping with IDH status and LGG status

 **Performance Evaluation:**

\- Confusion Matrix and Accuracy Metrics:

   \- High accuracy was achieved when separating IDH statuses.

   \- Some misclassifications occurred, particularly between closely related subtypes.

To further investigate the differences between IDH-mutant and IDH wild-type samples, a volcano plot was generated, highlighting up and downregulated CPG sites.
<img width="1465" height="677" alt="volcano" src="https://github.com/user-attachments/assets/fabf9db9-0fb3-4da3-968b-52bebeb5d333" />

Fig. 4: Volcano Plot for Mutant vs. Wild-Type

 **5\. Pathway Enrichment Analysis (Reactome)**

An overrepresentation analysis was performed using the Reactome database to identify pathways significantly enriched in genes associated with the CpG sites.

![pathway](https://github.com/user-attachments/assets/e7d4aa96-987e-47a6-8d63-d4f04cbae3af)


   Figure 5: Overrepresentation Analysis

**6\. Comparison with Target Paper**

The reference paper identified six distinct methylation clusters (LGm1-LGm6), separating IDH-mutant from wild-type tumors. Our analysis produced three clusters, which could be attributed to:

1\. Differences in CpG probe selection.

2\. Choice of clustering method (we used k-Means, while the paper used consensus clustering).

3\. Variations in preprocessing or dimensionality reduction parameters.

 **6\. Conclusion**

Our analysis using kNN demonstrated effective separation of glioma subtypes based on DNA methylation patterns. Although fewer clusters were identified compared to the target paper, the results align with known biological differences. Further work is required to achieve finer distinctions, potentially through improved probe selection and more advanced clustering techniques.

**7\. References**

\[1\]	N. Sharma *et al.*, “Isocitrate dehydrogenase mutations in gliomas: A review of current understanding and trials,” *Neuro-Oncology Adv.*, vol. 5, no. 1, Jan. 2023, doi: 10.1093/noajnl/vdad053.

\[2\]	V. Monga, K. Jones, and S. Chang, “CLINICAL RELEVANCE OF MOLECULAR MARKERS IN GLIOMAS,” *Rev. Médica Clínica Las Condes*, vol. 28, no. 3, pp. 343–351, May 2017, doi: 10.1016/j.rmclc.2017.05.003.

\[3\]	J. Fares, Y. Wan, R. Mair, and S. J. Price, “Molecular diversity in isocitrate dehydrogenase-wild-type glioblastoma,” *Brain Commun.*, vol. 6, no. 2, Mar. 2024, doi: 10.1093/braincomms/fcae108.
