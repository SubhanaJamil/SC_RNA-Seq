<h1 id="top"> 10X scRNA-seq Preprocessing Pipeline using Galaxy & STARsolo</h1>

<p>
A complete preprocessing workflow for 10X Genomics single-cell RNA-seq data using the Galaxy platform. 
This pipeline transforms raw FASTQ sequencing reads into a high-quality filtered count matrix 
ready for downstream single-cell analysis.
</p>

<ul>
  <li><a href="#overview">Overview</a></li>
  <li><a href="#tools">Tools Used</a></li>
  <li><a href="#dataset">Dataset</a></li>
  <li><a href="#workflow">Workflow Overview</a></li>

  <li><a href="#pipeline">Section: Preprocessing Pipeline</a>
    <ul>
      <li><a href="#upload">Data Upload</a></li>
      <li><a href="#chemistry">Identify Chemistry</a></li>
      <li><a href="#starsolo">Run STARsolo</a></li>
      <li><a href="#output">Inspect Output</a></li>
      <li><a href="#qc">Quality Control</a></li>
      <li><a href="#matrix">Raw Count Matrix</a></li>
      <li><a href="#filter">Cell Filtering</a></li>
      <li><a href="#rank">Barcode Rank Plot</a></li>
      <li><a href="#emptydrops">Advanced Filtering</a></li>
    </ul>
  </li>
</ul>

<h2 id="overview">Overview</h2>

<p>
This repository demonstrates preprocessing of 10X Genomics scRNA-seq data using Galaxy tools. 
The workflow follows best practices in single-cell transcriptomics, ensuring accurate 
demultiplexing, alignment, quantification, and filtering of cells before downstream analysis.
</p>

<ul>
  <li>Biological accuracy</li>
  <li>Efficient preprocessing</li>
  <li>Reproducibility</li>
  <li>Scalability</li>
</ul>

<h2 id="tools">Tools Used</h2>

<ul>
  <li><b>Galaxy:</b> Web-based platform for reproducible bioinformatics workflows</li>
  <li><b>STARsolo:</b> Alignment and quantification tool for scRNA-seq</li>
  <li><b>Cell Ranger (concept):</b> Reference pipeline for comparison</li>
  <li><b>MultiQC:</b> Quality control visualization</li>
  <li><b>DropletUtils:</b> Cell filtering and barcode analysis</li>
</ul>

<h2 id="dataset">Dataset</h2>

<ul>
  <li><b>Dataset:</b> PBMC 1K (subsampled)</li>
  <li><b>Platform:</b> 10X Genomics v3</li>
  <li><b>Input:</b> FASTQ files (R1 + R2)</li>
</ul>

<pre>
Files:
- R1 → Cell barcodes + UMI
- R2 → cDNA sequences
</pre>

<h2 id="workflow">Workflow Overview</h2>

<pre>
📥 FASTQ Data
      ↓
🧾 Upload & Organization
      ↓
🧬 STARsolo (Demultiplex + Align + Quantify)
      ↓
📊 Quality Control (MultiQC)
      ↓
📦 Raw Count Matrix
      ↓
🧹 Cell Filtering (DropletUtils)
      ↓
🧬 High-Quality Count Matrix
</pre>

<!-- PIPELINE SECTION -->

<h2 id="pipeline">🧪 Preprocessing Pipeline</h2>

<h3 id="upload">1️⃣ Data Upload</h3>

<p>
In this step, raw sequencing data and required reference files are uploaded into Galaxy. 
The FASTQ files contain sequencing reads where R1 stores cell barcode and UMI information, 
while R2 contains the actual gene expression sequences. Additionally, a gene annotation file 
and a barcode whitelist are required for accurate mapping and cell identification.
</p>

<pre><code>
Upload FASTQ files (R1, R2)
Upload GTF annotation file
Upload barcode whitelist
</code></pre>

<h3 id="chemistry">2️⃣ Identify 10X Chemistry</h3>

<p>
The chemistry version determines how barcode and UMI sequences are structured. 
This is identified by examining the length of reads in the R1 FASTQ file. 
Correct identification is essential because incorrect chemistry settings will 
lead to improper barcode extraction and inaccurate results.
</p>

<pre><code>
Check R1 read length

26 bp → v2 chemistry
28 bp → v3 chemistry
</code></pre>



<h3 id="starsolo">3️⃣ Run RNA STARsolo</h3>

<p>
STARsolo performs the core preprocessing steps including demultiplexing, alignment, 
and quantification. During this process, reads are assigned to cells using barcodes, 
aligned to the reference genome, and counted to produce gene expression values. 
UMI deduplication ensures that PCR duplicates are removed, improving accuracy.
</p>

<pre><code>
Reference Genome: hg19
Annotation: Homo_sapiens.GRCh37.75.gtf

Input:
- R1 → Barcode reads
- R2 → cDNA reads

Whitelist: 3M-february-2018.txt.gz
Chemistry: Chromium v3
UMI Deduplication: CellRanger-like
Filtering: Disabled
</code></pre>



<h3 id="output">4️⃣ Inspect Output Files</h3>

<p>
After processing, STARsolo generates multiple output files including alignment logs, 
BAM files, and count matrices. These outputs represent the raw gene expression data 
and mapping statistics, which must be evaluated before proceeding further.
</p>

<pre><code>
Output Files:
- Log file
- BAM file
- matrix.mtx
- barcodes.tsv
- genes.tsv
</code></pre>



<h3 id="qc">5️⃣ Quality Control (MultiQC)</h3>

<p>
Quality control is essential to assess sequencing and alignment performance. 
MultiQC aggregates STARsolo logs and provides visual summaries such as mapping rates 
and read quality. These metrics help determine whether the dataset is reliable.
</p>

<pre><code>
Tool: MultiQC
Input: STARsolo log file
</code></pre>



<h3 id="matrix">6️⃣ Raw Count Matrix</h3>

<p>
The raw count matrix contains gene expression values for all detected barcodes. 
However, this matrix includes many empty droplets and low-quality cells, making it 
highly sparse. Therefore, further filtering is required to obtain meaningful biological data.
</p>

<pre><code>
Matrix Components:
- matrix.mtx
- barcodes.tsv
- genes.tsv
</code></pre>


<h3 id="filter">7️⃣ Cell Filtering (DefaultDrops)</h3>

<p>
The DefaultDrops method removes low-quality cells and empty droplets based on statistical 
thresholds. This step reduces noise and produces a dataset similar to the output of the 
Cell Ranger pipeline.
</p>

<pre><code>
Tool: DropletUtils

Method: DefaultDrops
Expected Cells: 3000
Upper Quantile: 0.99
Lower Proportion: 0.1
</code></pre>


<h3 id="rank">8️⃣ Barcode Rank Plot</h3>

<p>
The barcode rank plot visualizes total UMI counts across all barcodes. 
The characteristic “knee point” separates real cells from empty droplets. 
This helps determine the number of high-quality cells in the dataset.
</p>

<pre><code>
Operation: Rank Barcodes
Lower Bound: 100
</code></pre>



<h3 id="emptydrops">9️⃣ Advanced Filtering (EmptyDrops)</h3>

<p>
EmptyDrops applies a statistical model to distinguish real cells from background noise. 
By controlling the false discovery rate (FDR), this method ensures that only high-confidence 
cells are retained for downstream analysis.
</p>

<pre><code>
Method: EmptyDrops
Lower Bound: 200
FDR: 0.01
</code></pre>



<h2>📊 Final Output</h2>

<p>
The final output is a filtered count matrix containing high-quality cells, 
typically around 270–300 cells. This dataset is now ready for downstream 
analysis such as clustering, visualization, and cell type annotation.
</p>

<h2>🚀 Next Steps</h2>

<p>
After preprocessing, the dataset can be analyzed using tools like Scanpy or Seurat. 
Typical downstream steps include normalization, dimensionality reduction (PCA, UMAP), 
clustering, and identification of marker genes.
</p>
