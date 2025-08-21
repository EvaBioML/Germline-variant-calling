# 🧬 ASD Germline Variant Calling Pipeline  

This project implements a **germline variant calling workflow** for an Autism Spectrum Disorder (ASD) study using NGS data.  
It covers raw data download, preprocessing, alignment, variant calling, filtering, and annotation.  

---
## 📂 Dataset  

- **NCBI SRA study:** [SRP187882](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP187882&o=acc_s%3Aa)  
- **Sample information:**  
  - `SRR8697627` – Mother  
  - `SRR8697637` – ASD child  
  - `SRR8697645` – Father
 

## 📂 Analysis Pipeline  

1. **Input Data**  
   - Raw sequencing reads: `FASTQ.gz` (paired-end)  

2. **Quality Control (QC)**  
   - Assess read quality with **FastQC**  
   - Adapter trimming and low-quality read removal (**Trimmomatic** )   

3. **Read Alignment**  
   - **DRAGEN** 
   - Convert SAM → BAM and sort using **SAMtools**  

4. **Post-alignment Processing**  
   - Mark duplicates (**Picard**)  
   - Base quality score recalibration (BQSR) with **GATK**  

5. **Variant Calling**  
   - Call germline variants using **GATK HaplotypeCaller**  

6. **Variant Filtering and Quality Control**  
   - Hard filtering   

7. **Annotation**  
   - Annotate variants with **ANNOVAR** ( **SnpEff**,  **VEP**  )
   - Include gene information, predicted effects, population frequencies (gnomAD), clinical significance (ClinVar)  

8. **Downstream Analysis**  
   - Identify rare / deleterious variants  
   - Family-based or cohort-based analysis  
   - Generate summary reports and variant statistics  

---

## ⚙️ Tools & Dependencies  

- **SRA Tools** (prefetch, grabseqs) – download and split SRA data  
- **FastQC** – quality control  
- **Trimmomatic** – adapter and quality trimming  
- **DRAGEN** – alignment and germline variant calling  
- **SAMtools** – BAM processing  
- **Picard** – sort BAM, mark duplicates  
- **GATK** – STR model calibration, HaplotypeCaller, VariantFiltration  
- **ANNOVAR** – variant annotation  
---




