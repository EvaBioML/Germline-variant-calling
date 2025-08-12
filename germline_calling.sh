#練習分析ASD，資料來源 https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP187882&o=acc_s%3Aa
#SRR8697627 mother  SRR8697637  autism  SRR8697645  father
# install bioconda::sra-tools
prefetch –O /PATH/ SRS4463419 SRS4463437 SRS4463427
#split
grabseqs sra -t 6 SRR8697627 SRR8697637 SRR8697645 

#Pre-trimming quality check 
# install bioconda::fastqc
fastqc -o /PATH/ -f fastq -t 48 /PATH/*.fastq.gz

#Trim adpater and low-quality positions  
#!/bin/bash
ADAPTERS=PATH/gatk4.2.6/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa
RAW_DIR=/PATH/
TRIM_DIR=/PATH/
TXT_FILE=/PATH//samples.txt

while read -r SAMPLE; do
    echo "正在處理樣本：$SAMPLE"

    R1="${RAW_DIR}/${SAMPLE}.fastq.gz"
    R2="${RAW_DIR}/${SAMPLE}_R2.fastq.gz"

    # 檢查檔案是否存在
    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        echo "⚠️  找不到 $SAMPLE 的 R1 或 R2 檔案，跳過..."
        continue
    fi

    trimmomatic PE -threads 12 \
        "$R1" "$R2" \
        "$TRIM_DIR/${SAMPLE}_P_R1.fastq.gz" \
        "$TRIM_DIR/${SAMPLE}_U_R1.fastq.gz" \
        "$TRIM_DIR/${SAMPLE}_P_R2.fastq.gz" \
        "$TRIM_DIR/${SAMPLE}_U_R2.fastq.gz" \
        ILLUMINACLIP:"$ADAPTERS":2:30:10:8:True \
        LEADING:20 TRAILING:20 MINLEN:36 \
        2> "$TRIM_DIR/${SAMPLE}_trimmomatic.log"

done echo "完成"


#Part 2. Germline variant calling pipeline in DRAGEN mode
#dragen-os --build-hash-table true --ht-reference /PATH/GRCh38.d1.vd1.fa  --output-directory /PATH/dragmap/
#gatk ComposeSTRTableFile -R /PATH/GRCh38.d1.vd1.fa -O /PATH/GRCh38.d1.vd1.str.table.tsv
# Alignment (Input:trimmed fastq file, Output: raw bam file, Tool: dramap)(做完3人)
dragen-os \
-r /PATH/hg38/ \
--RGID ASD \
--RGSM mom \
-1 /PATH/SRR8697627_1.fastq \
-2 /PATH/SRR8697627_2.fastq | \
samtools view \
-@ 48 \
-b -h \
-o /PATH/mom.bam


# Sort bam and mark duplicates (Input:raw bam file, Output: sorted and marked bam file, Tool: picard)(做完3人) 

picard "-Xmx3G -Djava.io.tmpdir=./"  SortSam \
CREATE_INDEX=true \
INPUT=/your path/asd.bam \
OUTPUT=/your path/asd_sorted.bam \
SORT_ORDER=coordinate \
VALIDATION_STRINGENCY=STRICT

picard "-Xmx3G -Djava.io.tmpdir=./"  MarkDuplicates \
CREATE_INDEX=true \
INPUT=/your path/asd_sorted.bam \
OUTPUT=/your path/asd_sorted_markdup.bam \
METRICS_FILE=/your path/asd_sorted_markdup_metrics.txt \
VALIDATION_STRINGENCY=STRICT

# Calibrate DRAGEN STR model  (Input:sorted and marked bam file, Output: STR recalibration model file, Tool: GATK4) 
# Germline variant calling and hard filtering (variant filtering) on "QUAL"  (Input:sorted and marked bam file, Output: filterd germline VCF, Tool: GATK4) 
mkdir -p /your path/BamOut/
mkdir -p /your path/Analysis_Reday/

Gatk CalibrateDragstrModel \
-R path/hg38/Homo_sapiens_assembly38.fasta \
-I path/asd_sorted_markdup.bam \
-str path/dragmap_hg38/Homo_sapiens_assembly38.str.table.tsv \
-O path/asd_dragstr_model.txt


gatk HaplotypeCaller \
-R path/ref/GRCh38.d1.vd1.fa \
-I path/asd_sorted_markdup.bam \
-L path/asd_Covered.bed \
-O path/asd.germline.vcf \
--dragen-mode true \
--dragstr-params-path path/asd_dragstr_model.txt \
-bamout path/BamOut/asd.bamout.bam

gatk VariantFiltration \
-V path/asd.germline.vcf \
--filter-expression "QUAL < 10.4139" \
--filter-name "DRAGENHardQUAL" \
-O path/asd.germline.hf.vcf

gatk SelectVariants \
--java-options "-Xmx5G" \
-R path/GRCh38.d1.vd1.fa \
--exclude-filtered \
-L path/asd_Covered.bed \
-V  path/asd.germline.hf.vcf \
-O  path/asd.germline.final.vcf


##annotation
python Intervar.py -b hg38 -t path/intervardb/ -i path/prefix_asd.avinput  --table_annovar=path/annovar/table_annovar.pl --convert2annovar=path/annovar/convert2annovar.pl --annotate_variation=path/annovar/annotate_variation.pl --database_locat=path//humandb/ --input_type=AVinput  -o path/annovar/prefix_asd
/your annovar path/table_annovar.pl /your VCF path/your.vcf \
path/humandb/ -buildver hg38 -out /your output path/prefix_name -remove -protocol \
refGene,\
ensGene,\
cytoBand,\
wgRna,\
genomicSuperDups,\
rmsk,\
avsnp150,\
cosmic96_coding,\
cosmic96_noncoding,\
clinvar_20220320,\
1000g2015aug_eas,\
1000g2015aug_sas,\
1000g2015aug_eur,\
1000g2015aug_afr,\
1000g2015aug_amr,\
1000g2015aug_all,\
gnomad211_exome,\
esp6500siv2_all,\
exac03nontcga,\
TaiwanBiobank993WGS,\
dbnsfp42a,\
dbscsnv11,\
gwasCatalog \
-operation \
g,g,r,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r -nastring . -vcfinput


python /your InterVar path/Intervar.py \
-b hg38 \
-t path/intervardb/ \
-i /your variant avinput format path/your.avinput  \
--table_annovar=/your annovar path/table_annovar.pl \
--convert2annovar=/your annovar path/convert2annovar.pl \
--annotate_variation=/your annovar path/annotate_variation.pl \
--database_locat= path/humandb/ \
--input_type=AVinput  \
-o /your output path/prefix_name
