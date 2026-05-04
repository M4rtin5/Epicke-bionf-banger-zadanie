# Analýza sekvenačných dát – S14 (tumor vs control)

## 1. Kontrola surových dát

Pomocou `seqkit stats` boli analyzované FASTQ súbory:

* **Tumor (S14.T)**:

  * počet čítaní: 60,346,220
  * priemerná dĺžka: ~101 bp

* **Control (S14.C)**:

  * počet čítaní: 46,678,977
  * priemerná dĺžka: ~101 bp

Dáta sú konzistentné a zodpovedajú Illumina sekvenovaniu.

Kvalita bola overená nástrojom **FastQC**.

Commandy:

```
fastqc S14.T_R1.fastq.gz S14.T_R2.fastq.gz -o fastqc_bombo/
fastqc S14.C_R1.fastq.gz S14.C_R2.fastq.gz -o fastqc_nezamestnany/
```

---

## 2. Mapovanie na referenčný genóm

Mapovanie bolo vykonané pomocou `bwa mem` na referenciu hg38 a následne zoradené pomocou `samtools sort`.
Celý príkaz:

```
bwa mem -t 2 hg38.fa S14.T_R1.fastq.gz S14.T_R2.fastq.gz | samtools sort -@ 2 -m 512M -o S14.T.sorted.bam
bwa mem -t 2 hg38.fa S14.C_R1.fastq.gz S14.C_R2.fastq.gz | samtools sort -@ 2 -m 512M -o S14.C.sorted.bam
```

Indexovanie:

```
samtools index S14.T.sorted.bam
samtools index S14.C.sorted.bam
```

---

## 3. Pridanie read groups

Pomocou GATK (`AddOrReplaceReadGroups`) boli pridané read group informácie.

```
./gatk-4.4.0.0/gatk AddOrReplaceReadGroups -I S14.T.sorted.bam -O S14.T.rg.bam -RGID S14T -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM TUMOR
./gatk-4.4.0.0/gatk AddOrReplaceReadGroups -I S14.C.sorted.bam -O S14.C.rg.bam -RGID S14C -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM NORMAL
```

---

## 4. Kontrola mapovania

Vypísanie štatistík koľko percent čítaní sa podarilo prilepiť na genóm.

```
samtools index S14.T.rg.bam
samtools index S14.C.rg.bam
```

Použitý nástroj `samtools flagstat`:

* **Tumor:**

  * mapped: 99.78 %
  * properly paired: 99.41 %

* **Control:**

  * mapped: 99.80 %
  * properly paired: 99.42 %

Mapovanie bolo veľmi úspešné.

```
samtools flagstat S14.T.rg.bam > S14.T.rg.flagstat.txt
samtools flagstat S14.C.rg.bam > S14.C.rg.flagstat.txt
```

---

## 5. Coverage a kvalita

* `samtools coverage S14.T.rg.bam > S14.T.coverage.txt` – výstup bol rozsiahly
* `samtools coverage S14.C.rg.bam > S14.C.coverage.txt` – výstup bol rozsiahly
* `qualimap bamqc -bam S14.C.rg.bam -outdir qualimap_report_C` a `qualimap bamqc -bam S14.T.rg.bam -outdir qualimap_report_T --java-mem-size=6G -nt 4` – vizualizácia:

  * coverage
  * mapping quality
  * distribúcia fragmentov

---

## 6. Odstránenie duplicít

Použitý GATK nástroj:

```
./gatk-4.4.0.0/gatk MarkDuplicates -I S14.T.rg.bam -O S14.T.rg.dedup.bam -M T.dup_metrics.txt --CREATE_INDEX
./gatk-4.4.0.0/gatk MarkDuplicates -I S14.C.rg.bam -O S14.C.rg.dedup.bam -M C.dup_metrics.txt --CREATE_INDEX
```

Výstup:

* deduplikované BAM súbory
* metriky v `.txt` súboroch

---

## 7. Base Quality Score Recalibration (BQSR)

Použitá databáza dbSNP.

Kroky:

1. `BaseRecalibrator`
2. `ApplyBQSR`

Krok 1:
```
./gatk-4.4.0.0/gatk BaseRecalibrator -R hg38.fa -I S14.T.rg.dedup.bam --known-sites dbsnp.chr.vcf.gz -O T.recal_data.recalib
./gatk-4.4.0.0/gatk BaseRecalibrator -R hg38.fa -I S14.C.rg.dedup.bam --known-sites dbsnp.chr.vcf.gz -O C.recal_data.recalib
```

Krok 2:
```
./gatk-4.4.0.0/gatk ApplyBQSR -R hg38.fa -I S14.T.rg.dedup.bam --bqsr-recal-file T.recal_data.recalib -O T.bwa.cleaned.bam
./gatk-4.4.0.0/gatk ApplyBQSR -R hg38.fa -I S14.C.rg.dedup.bam --bqsr-recal-file C.recal_data.recalib -O C.bwa.cleaned.bam
```

Výsledkom sú recalibrované BAM súbory:

* `T.bwa.cleaned.bam`
* `C.bwa.cleaned.bam`

---

## 8. Variant calling

Použitý nástroj:
* `GATK Mutect2`

Commandy na variant calling:

```
NORMAL_SAMPLE=$(samtools samples C.bwa.cleaned.bam | cut -f1)
./gatk-4.4.0.0/gatk Mutect2 -R hg38.fa -I T.bwa.cleaned.bam -I C.bwa.cleaned.bam -normal $NORMAL_SAMPLE -O Tumor_Control.vcf.gz
```

Vstupy:

* tumor + control BAM
* referenčný genóm hg38

Výstup:

* `Tumor_Control.vcf.gz`

---

## 9. Filtrovanie variantov

Použitý nástroj:
* `GATK Mutect2`

Command na filtrovanie:
```
./gatk-4.4.0.0/gatk FilterMutectCalls -R hg38.fa -V Tumor_Control.vcf.gz -O Tumor_Control.filtered.vcf.gz
```

Výsledok:

* `Tumor_Control.filtered.vcf.gz`

---

## 10. Detekcia štrukturálnych variantov (Manta)

Na identifikáciu štrukturálnych variantov (SV) bol použitý nástroj Manta.

Konfigurácia:

```bash
configManta.py \
  --tumorBam T.bwa.cleaned.bam \
  --normalBam C.bwa.cleaned.bam \
  --referenceFasta hg38.fa \
  --runDir manta_run
```

Spustenie analýzy:
```bash
manta_run/runWorkflow.py -m local -j 4
```

Hlavný výsledný súbor:

* `results/variants/somaticSV.vcf.gz`

Obsahuje:

* delécie (DEL)
* inzercie (INS)
* duplikácie (DUP)
* inverzie (INV)

---

## Záver

Pipeline úspešne prešiel všetkými krokmi:

kvalitné mapovanie (~99.8 %)
správne spracovanie BAM súborov
úspešné variant calling a filtrácia

Výsledkom je finálny VCF súbor obsahujúci somatické mutácie medzi tumor a control vzorkou.
