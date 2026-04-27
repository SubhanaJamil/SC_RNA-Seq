<h1> <p align="center"> Single Cell RNA-Seq </p></h1>
<h2>📌 Table of Contents</h2>

<ul>
  <li><a href="#overview">Overview</a></li>
  <li><a href="#tools-used">Tools Used</a></li>
  <li><a href="#dataset">Dataset</a></li>
  <li><a href="#workflow-overview">Workflow Overview</a></li>
  <li><a href="#installation">Installation</a></li>

  <li><b>Pipeline Implementation</b>
    <ul>
      <li><a href="#data-loading">Data Loading</a></li>
      <li><a href="#quality-control">Quality Control</a></li>
      <li><a href="#quality-visualization">Quality Visualization</a></li>
      <li><a href="#filtering">Filtering</a></li>
      <li><a href="#doublet-detection">Doublet Detection</a></li>
      <li><a href="#normalization">Normalization</a></li>
      <li><a href="#hvg-selection">Highly Variable Genes</a></li>
      <li><a href="#pca">PCA</a></li>
      <li><a href="#neighbors-umap">Neighbors & UMAP</a></li>
      <li><a href="#clustering">Clustering</a></li>
      <li><a href="#qc-recheck">QC Recheck</a></li>
      <li><a href="#marker-genes">Marker Genes</a></li>
      <li><a href="#dotplot">Dotplot Visualization</a></li>
      <li><a href="#cell-annotation">Cell Annotation</a></li>
      <li><a href="#differential-expression">Differential Expression</a></li>
    </ul>
  </li>
  
  <li><a href="#results">Results (Top Differentially Expressed Genes)</a></li>
  <li><a href="#interpretation">Results Interpretation</a></li>
  <li><a href="#conclusion">Conclusion</a></li>
</ul>
<h2 id="overview">Overview</h2>
<p>
This project implements a complete scRNA-seq preprocessing, clustering, and cell-type annotation pipeline using publicly available bone marrow mononuclear cell data. The workflow follows best practices in single-cell transcriptomics, ensuring biological accuracy, reproducibility, and scalability. It integrates statistical rigor with computational efficiency to extract meaningful cellular heterogeneity.
</p>

<p>
The workflow follows <b>best practices in single-cell transcriptomics</b>, ensuring:
</p>

<ul>
  <li>Biological accuracy</li>
  <li>Statistical robustness</li>
  <li>Reproducibility</li>
  <li>Scalability</li>
</ul>

<h2>Tools Used</h2>

<p>
This pipeline was implemented using a standard single-cell RNA-seq computational ecosystem, primarily built on Scanpy + AnnData, with Pooch for data retrieval, Scrublet for doublet detection, and Matplotlib for visualization.:
</p>

<ul>
  <li><b>Scanpy: </b>  End-to-end single-cell analysis framework for preprocessing, dimensionality reduction, clustering, and visualization.</li>
  <li><b>AnnData: </b> Efficient data structure for storing gene expression matrices with annotations.</li>
  <li><b>Pooch: </b>  Lightweight data retrieval system for reproducible dataset downloading and caching.</li>
  <li><b>NumPy: </b>  High-performance numerical computations and array manipulation.</li>
  <li><b>Matplotlib: </b> Core visualization library for generating publication-quality plots.</li>
  <li><b>Scrublet: </b> Algorithm for identifying doublets in droplet-based single-cell sequencing data.</li>
</ul>

<h2 id="dataset">Dataset</h2>

<ul>
  <li><b>Source:</b> NeurIPS 2021 Benchmark Dataset</li>
  <li><b>Platform:</b> 10X Multiome (Gene Expression + Chromatin Accessibility)</li>
  <li><b>Samples used:</b> s1d1, s1d3</li>
</ul>

<p><b>Final Dataset Dimensions:</b></p>

<pre>
AnnData object with n_obs × n_vars = 17125 × 36601
</pre>

<h2 id="workflow-overview">Workflow Overview</h2>

<pre>

╭──── ──────────────── 🧬 scRNA-seq PIPELINE ───────────────────────╮
                      
                       📥 Raw Data Acquisition
                                   ↓
                       🧾 Data Loading & Merging
                                   ↓
                      🧪 Quality Control (QC Metrics + Visualization)
                                   ↓
                      🧹 Filtering (Cells & Genes)
                                   ↓
                      🧬 Doublet Detection (Scrublet)
                                   ↓
                      📏 Normalization & Log Transformation
                                   ↓
                      ⭐ Highly Variable Gene Selection
                                   ↓
                      📉 PCA (Dimensionality Reduction)
                                   ↓
                      🌐 Graph Construction (Neighbors)
                                   ↓
                      🗺️ UMAP Visualization
                                   ↓
                      🧩 Clustering (Leiden Algorithm)
                                   ↓
                      🧬 Marker Gene Identification
                                   ↓
                      🏷️ Cell Type Annotation
                                   ↓
                      📊 Differential Expression Analysis
   

╰───────────────────── 🚀 END OF PIPELINE ───────────────────────────╯
</pre>

---

<h2 id="installation">Installation</h2>

<p>Install the required dependencies using pip:</p>

<pre>
pip install scanpy anndata pooch
</pre>

<h2 id="pipeline-implementation">Pipeline Implementation</h2>

<h3 id="data-loading">1️⃣ Data Loading</h3>

<p>
This step retrieves raw sequencing data using pooch, loads it into AnnData objects, and merges multiple samples into a single dataset. Ensuring unique gene and cell identifiers is critical to avoid downstream conflicts.
</p>

<pre><code class="language-python">
import scanpy as sc
import anndata as ad
import pooch

sc.set_figure_params(dpi=80, facecolor="white")

DATA = pooch.create(
    path=pooch.os_cache("scRNA_data"),
    base_url="doi:10.6084/m9.figshare.22716739.v1/",
)
DATA.load_registry_from_doi()

samples = {
    "s1d1": "s1d1_filtered_feature_bc_matrix.h5",
    "s1d3": "s1d3_filtered_feature_bc_matrix.h5",
}

adatas = {}
for sid, file in samples.items():
    path = DATA.fetch(file)
    adata = sc.read_10x_h5(path)
    adata.var_names_make_unique()
    adatas[sid] = adata

adata = ad.concat(adatas, label="sample")
adata.obs_names_make_unique()
adata.var_names_make_unique()
print(adata)

</code></pre>

<h3 id="quality-control">2️⃣ Quality Control</h3>

<p>
Quality control quantifies technical artifacts such as mitochondrial contamination and sequencing depth. These metrics help distinguish high-quality cells from stressed or dying cells. It Computes QC metrics to identify low-quality or stressed cells.
</p>

<pre><code class="language-python">
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
adata.var["hb"] = adata.var_names.str.contains("^HB")

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt", "ribo", "hb"],
    inplace=True,
    log1p=True
)
</code></pre>

<h3 id="quality-visualization">3️⃣ Quality Visualization</h3>
<p> Visualization allows intuitive identification of outliers and abnormal cells. Patterns such as high mitochondrial percentage or low gene counts indicate poor-quality cells. </p>

<h4>Scatter Plot</h4>
<p>
This plot shows key quality metrics across all cells including gene counts, total counts, and mitochondrial gene percentage.
It helps identify low-quality or stressed cells.
</p>

<pre><code class="language-python">
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    multi_panel=True
)
</code></pre>

<img width="2269" height="950" alt="Figure_1" src="https://github.com/user-attachments/assets/83bdfb6b-a636-4779-ad65-8c2c3ac17b50" />

<h4>Scatter Plot</h4>

<p>
This scatter plot visualizes the relationship between total counts and number of detected genes, colored by mitochondrial content.
It helps detect outliers and poor-quality cells.
</p>

<pre><code class="language-python">
sc.pl.scatter(
    adata,
    x="total_counts",
    y="n_genes_by_counts",
    color="pct_counts_mt"
)
</code></pre>
<p align="center">
  <img width="454" height="400" alt="Figure_2" src="https://github.com/user-attachments/assets/657a158e-fe74-4e9d-9fe4-320d9f9c3fe6" />
</p>



<h3 id="filtering">4️⃣ Filtering</h3>
<p>
Filtering removes noise by excluding low-quality cells and rarely expressed genes. This step improves downstream statistical power, clustering & analysis accuracy.
</p>

<pre><code class="language-python">
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)
</code></pre>



<h3 id="doublet-detection">5️⃣ Doublet Detection</h3>
<p> Doublets arise when two cells are captured as one, leading to misleading hybrid expression profiles. Detecting and flagging them prevents incorrect biological interpretation.</p>

<pre><code class="language-python">
sc.pp.scrublet(adata, batch_key="sample")
</code></pre>


<h3 id="normalization">6️⃣ Normalization</h3>
<p> Normalization corrects for sequencing depth differences across cells, while log transformation stabilizes variance and improves comparability across genes. </p>

<pre><code class="language-python">
adata.raw = adata.copy()
adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
</code></pre>

<h3 id="hvg-selection">7️⃣ Highly Variable Genes</h3>
<p> Only the most informative genes are retained to reduce dimensionality and noise. HVGs capture biological variability critical for clustering. </p>
<pre><code class="language-python">
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    batch_key="sample"
)
sc.pl.highly_variable_genes(adata)
adata = adata[:, adata.var.highly_variable].copy()
</code></pre>

<img width="1200" height="600" alt="Figure_3" src="https://github.com/user-attachments/assets/96170d50-6e1c-4f3f-9f75-a090f3c6e934" />


<h3 id="pca">8️⃣ Dimensionality Reduction (PCA)</h3>
<P> Principal Component Analysis compresses high-dimensional data into fewer dimensions while preserving variance. This step forms the basis for clustering. </P>

<pre><code class="language-python">
sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50)
  
</code></pre>
<p align="center">
<img width="445" height="400" alt="Figure_4" src="https://github.com/user-attachments/assets/f8962303-5112-451f-88a2-b4fa7b2862a6" />
</p>

<h3 id="neighbors-umap">9️⃣ Graph Construction & UMAP</h3>
<p> A nearest-neighbor graph models cell similarity, and UMAP projects this structure into a 2D space for visualization of cellular relationships. </p>

<pre><code class="language-python">
sc.pp.neighbors(adata, n_pcs=30)
sc.tl.umap(adata)
  
sc.pl.umap(adata, color="sample")
</code></pre>
<p align="center">
<img width="445" height="400" alt="Figure_5" src="https://github.com/user-attachments/assets/317fd8c3-3846-4852-b8c9-8e91c49d79e4" />
</p>

<h3 id="clustering">🔟 Clustering (Leiden)</h3>
<p> Clustering groups cells into biologically meaningful populations. Multiple resolutions provide hierarchical insights from broad to fine-grained clusters. </p>

<pre><code class="language-python">
sc.tl.leiden(adata, resolution=0.5, flavor="igraph")

for res in [0.02, 0.5, 2.0]:
    sc.tl.leiden(adata, key_added=f"leiden_res_{res}", resolution=res, flavor="igraph")

sc.pl.umap(
    adata,
    color=["leiden_res_0.02", "leiden_res_0.5", "leiden_res_2.0"]
)
</code></pre>
<p align="center">
<img width="445" height="400" alt="Figure_6" src="https://github.com/user-attachments/assets/e8210ef7-adaa-4bfb-9e82-1d2b48e87137" />
</p>

<h3 id="qc-recheck">1️⃣1️⃣ QC Recheck</h3>
<p> Post-clustering QC ensures clusters are biologically meaningful and not driven by artifacts such as doublets or mitochondrial bias. </p>
<pre><code class="language-python">
sc.pl.umap(
    adata,
    color=["leiden", "doublet_score", "pct_counts_mt"]
)
</code></pre>

<img width="2172" height="600" alt="Figure_7" src="https://github.com/user-attachments/assets/47cff3bb-8a76-4c87-b8e2-a0ee9d7ebbb1" />


<h3 id="marker-genes">1️⃣2️⃣ Marker Genes</h3>
<p> Marker genes define cell identity. These canonical markers are used to interpret clusters biologically. </p>
<pre><code class="language-python">
marker_genes = {
    "CD4+ T": ["CD4", "IL7R"],
    "CD8+ T": ["CD8A", "CD8B"],
    "NK": ["GNLY", "NKG7"],
    "B cells": ["MS4A1"],
    "Monocytes": ["CD14", "LYZ"]
}
</code></pre>
<img width="995" height="598" alt="Figure_8(b)" src="https://github.com/user-attachments/assets/071ca897-6452-4059-b0b4-bf4e2c4d9ec3" />

<h3 id="dotplot">1️⃣3️⃣ Dotplot Visualization</h3>
<p> Dotplots summarize gene expression across clusters, combining expression level and frequency into a single intuitive visualization. </p>
<pre><code class="language-python">
sc.tl.dendrogram(adata, groupby="leiden_res_0.5")

sc.pl.dotplot(
    adata,
    marker_genes,
    groupby="leiden_res_0.5",
    use_raw=True,
    standard_scale="var"
)
</code></pre>

<h3 id="cell-annotation">1️⃣4️⃣ Cell Type Annotation</h3>
<p> Clusters are mapped to known biological cell types using marker gene expression patterns, translating computational clusters into biological meaning. </p>

<pre><code class="language-python">
adata.obs["cell_type"] = adata.obs["leiden_res_0.5"].map({
    "0": "T cells",
    "1": "Monocytes",
    "2": "B cells",
    "3": "NK cells"
})
</code></pre>

<img width="995" height="587" alt="Figure_8" src="https://github.com/user-attachments/assets/d56b2d53-814e-4bd7-84e8-c0f04a3e0af5" />

<h3 id="differential-expression">1️⃣5️⃣ Differential Expression</h3>
<p> Differential expression identifies genes that distinguish clusters, revealing functional differences and validating cell identities.</p>

<pre><code class="language-python">
sc.tl.rank_genes_groups(
    adata,
    groupby="leiden_res_0.5",
    method="wilcoxon"
)
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5)

print(sc.get.rank_genes_groups_df(adata, group="0").head(5))
</code></pre>

<img width="2885" height="1095" alt="Figure_9" src="https://github.com/user-attachments/assets/7c47a7e5-7297-401b-99b7-b4adc213a585" />

<h2 id="results">Results (Top Differentially Expressed Genes)</h2>
<p>
The following genes were identified as top differentially expressed markers using Wilcoxon rank-sum testing. These genes strongly define major immune cell populations, especially T cells.
</p>

<pre>
┌─────────┬────────┬────────────────┬────────┬────────────┐
│ Gene    │ Score  │ LogFC          │ p-val  │ adj p-val  │
├─────────┼────────┼────────────────┼────────┼────────────┤
│ CD3D    │ 77.57  │ 6.50           │ 0.0    │ 0.0        │
│ CD3E    │ 74.02  │ 5.18           │ 0.0    │ 0.0        │
│ BCL11B  │ 73.44  │ 5.80           │ 0.0    │ 0.0        │
│ CD3G    │ 73.06  │ 5.11           │ 0.0    │ 0.0        │
│ RPL30   │ 72.25  │ 54.56          │ 0.0    │ 0.0        │
└─────────┴────────┴────────────────┴────────┴────────────┘
</pre>

<h2 id="interpretation">Results Interpretation</h2>
<p>
Top marker genes such as <b>CD3D, CD3E, and CD3G</b> confirm accurate identification of <b>T-cell populations</b>.
</p>

<p>
Extremely low p-values (~0) reflect <b>high statistical power</b> due to large sample size, not computational error.
</p>


<h2 id="conclusion">Conclusion</h2>

<p>
This pipeline successfully identifies major immune cell populations:
</p>

<ul>
  <li>🧬 T cells</li>
  <li>🧬 B cells</li>
  <li>🧬 NK cells</li>
  <li>🧬 Monocytes</li>
</ul>

<p>
It provides a <b>robust, scalable foundation</b> for:
</p>

<ul>
  <li>Disease modeling</li>
  <li>Biomarker discovery</li>
  <li>Multi-omics integration</li>
</ul>

---
