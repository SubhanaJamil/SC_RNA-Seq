

<h1> Single-Cell RNA-seq Preprocessing Pipeline (STARsolo)</h1>

<section>
<h2>1️⃣ Quality Control (FastQC + MultiQC)</h2>
<pre>
mkdir -p qc_reports

fastqc R1.fastq.gz R2.fastq.gz -o qc_reports/

multiqc qc_reports/ -o qc_reports/
</pre>
<p>This step evaluates raw sequencing quality using FastQC and summarizes results using MultiQC.</p>
</section>

<section>
<h2>2️⃣ STAR Genome Index Generation (hg19)</h2>
<pre>
STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir hg19_index \
--genomeFastaFiles hg19.fa \
--sjdbGTFfile annotation.gtf \
--sjdbOverhang 100
</pre>
<p>This step builds the reference genome index required for alignment.</p>
</section>

<section>
<h2>3️⃣ STARsolo Alignment + Quantification</h2>
<pre>
STAR --runThreadN 8 \
--genomeDir hg19_index \
--readFilesIn R2.fastq.gz R1.fastq.gz \
--readFilesCommand zcat \
--outFileNamePrefix output/star_solo_ \
--outSAMtype BAM SortedByCoordinate \
--soloType CB_UMI_Simple \
--soloCBlen 16 \
--soloUMIlen 12 \
--soloCBwhitelist whitelist.txt \
--soloFeatures Gene
</pre>
<p>This step performs alignment and generates gene-cell count matrices using STARsolo.</p>
</section>

<section>
<h2>4️⃣ Output Structure (STARsolo)</h2>
<pre>
output/star_solo_Solo.out/Gene/filtered/

Key files:
- matrix.mtx
- barcodes.tsv
- features.tsv
</pre>
<p>Filtered expression matrix is generated for downstream analysis.</p>
</section>

<section>
<h2>5️⃣ MultiQC Summary Report</h2>
<pre>
multiqc output/ -o multiqc_report/
</pre>
<p>Aggregates QC reports into a single interactive HTML report.</p>
</section>

<section>
<h2>6️⃣ R: Load Matrix + QC Analysis</h2>
<pre>
library(DropletUtils)
library(Matrix)
library(scater)

sce <- read10xCounts("output/star_solo_Solo.out/Gene/filtered/")
counts_matrix <- counts(sce)

sce$mito_percent <- calculateQCMetrics(sce)$subsets_Mito_percent
</pre>
<p>Loads STARsolo output and computes QC metrics.</p>
</section>

<section>
<h2>7️⃣ Barcode Rank Plot (Knee Plot)</h2>
<pre>
br <- barcodeRanks(counts_matrix)

plot(br$rank, br$total,
     log="xy",
     xlab="Rank",
     ylab="Total UMI Counts",
     main="Barcode Rank Plot")

abline(v=metadata(br)$knee, col="blue", lty=2)
abline(v=metadata(br)$inflection, col="red", lty=2)
</pre>
<p>Identifies true cells vs empty droplets.</p>
</section>

<section>
<h2>8️⃣ Cell Filtering (EmptyDrops)</h2>
<pre>
set.seed(123)

ed <- emptyDrops(counts_matrix)

filtered_cells <- which(ed$FDR <= 0.01)

filtered_matrix <- counts_matrix[, filtered_cells]

length(filtered_cells)
</pre>
<p>Removes empty droplets and retains high-confidence cells.</p>
</section>

<section>
<h2>9️⃣ Basic Quality Filtering (Optional)</h2>
<pre>
keep <- sce$total_counts > 500 &
        sce$total_features_by_counts > 200 &
        sce$mito_percent < 10

sce_filtered <- sce[, keep]
</pre>
<p>Applies thresholds to remove low-quality cells.</p>
</section>

