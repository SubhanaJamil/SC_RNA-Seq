<h1 id="top">🧬 Getting Started with AnnData & Scanpy (Single-Cell Data Structure Tutorial)</h1>

<p>
A complete beginner-friendly guide to understanding the <b>AnnData object structure</b> used in Scanpy and scverse for single-cell RNA-seq analysis.
This tutorial explains how single-cell data is stored, accessed, and manipulated in Python.
</p>

<ul>
  <li><a href="#overview">Overview</a></li>
  <li><a href="#tools-used">Tools Used</a></li>
  <li><a href="#workflow">Workflow Overview</a></li>
  <li><a href="#anndata-object">AnnData Object Structure</a>
    <ul>
      <li><a href="#x-matrix">Active Data Matrix (X)</a></li>
      <li><a href="#layers">Layers</a></li>
      <li><a href="#obs-var">Cell & Gene Annotations (obs & var)</a></li>
      <li><a href="#indexing">Indexing (obs_names & var_names)</a></li>
      <li><a href="#obsm">Multidimensional Data (obsm)</a></li>
      <li><a href="#obsp">Pairwise Data (obsp)</a></li>
      <li><a href="#uns">Unstructured Metadata (uns)</a></li>
    </ul>
  </li>
  <li><a href="#subsetting">Subsetting AnnData</a></li>
  <li><a href="#views-copies">Views vs Copies</a></li>
  <li><a href="#conclusion">Conclusion</a></li>
</ul>

---

<h2 id="overview">📌 Overview</h2>

<p>
This tutorial demonstrates how single-cell RNA-seq data is stored in the <b>AnnData</b> object, the core data structure used in Scanpy.
It helps users understand how to explore datasets, metadata, and embeddings efficiently.
</p>

<ul>
  <li>Understanding AnnData structure</li>
  <li>Accessing expression matrices</li>
  <li>Working with annotations</li>
  <li>Handling embeddings and graphs</li>
</ul>

---

<h2 id="tools-used">🧪 Tools Used</h2>

<ul>
  <li><b>anndata</b> – Core data structure for single-cell analysis</li>
  <li><b>scanpy</b> – Single-cell analysis toolkit</li>
  <li><b>numpy</b> – Numerical computations</li>
  <li><b>matplotlib</b> – Visualization</li>
  <li><b>pooch</b> – Data downloading and caching</li>
</ul>

---

<h2 id="workflow">⚙️ Workflow Overview</h2>

<pre>
Raw .h5ad file → Load with anndata → Explore structure → Analyze components → Subset data → Handle views/copies
</pre>

---

<h2 id="anndata-object">🧬 AnnData Object Structure</h2>

<p>
AnnData stores single-cell data in a structured format combining expression data and metadata.
</p>

---

<h3 id="x-matrix">📊 1. Active Data Matrix (X)</h3>

<p>
The <code>adata.X</code> matrix contains the main gene expression data (cells × genes).
It is usually stored as a sparse matrix for memory efficiency.
</p>

---

<h3 id="layers">📦 2. Layers</h3>

<p>
Layers store alternative versions of the expression matrix (e.g., raw counts, normalized counts).
</p>

<pre>
adata.layers["raw"]
adata.layers["normalized"]
</pre>

---

<h3 id="obs-var">🧫 3. Cell & Gene Annotations (obs & var)</h3>

<ul>
  <li><b>obs</b> → cell-level metadata (QC metrics, clusters)</li>
  <li><b>var</b> → gene-level metadata (gene names, IDs)</li>
</ul>

---

<h3 id="indexing">🔑 4. Indexing (obs_names & var_names)</h3>

<p>
Each cell and gene has a unique identifier for consistent referencing across analysis.
</p>

<pre>
adata.obs_names  # cell IDs
adata.var_names  # gene names or IDs
</pre>

---

<h3 id="obsm">📉 5. Multidimensional Data (obsm)</h3>

<p>
Stores embeddings like PCA, t-SNE, and UMAP.
</p>

<pre>
adata.obsm["X_pca"]
adata.obsm["X_umap"]
adata.obsm["X_tsne"]
</pre>

---

<h3 id="obsp">🔗 6. Pairwise Data (obsp)</h3>

<p>
Stores cell-cell similarity or distance matrices (e.g., nearest neighbor graphs).
</p>

<pre>
adata.obsp["distances"]
</pre>

---

<h3 id="uns">🗂️ 7. Unstructured Metadata (uns)</h3>

<p>
Stores dataset-level metadata such as clustering parameters, colors, and PCA variance.
</p>

<pre>
adata.uns["pca"]
adata.uns["louvain"]
</pre>

---

<h2 id="subsetting">✂️ Subsetting AnnData</h2>

<p>
AnnData supports flexible subsetting by:
</p>

<ul>
  <li>Index names (genes/cells)</li>
  <li>Numerical indexing</li>
  <li>Boolean filtering</li>
</ul>

<pre>
adata_subset = adata[:10, ["CD3D", "LYZ"]]
</pre>

---

<h2 id="views-copies">🔄 Views vs Copies</h2>

<p>
Subsetting creates a <b>view</b> (not a copy), meaning it references original data.
</p>

<ul>
  <li>View → memory efficient</li>
  <li>Copy → independent dataset</li>
</ul>

<pre>
adata_view = adata[:5, :]
adata_copy = adata_view.copy()
</pre>

---

<h2 id="conclusion">🏁 Conclusion</h2>

<p>
AnnData is the core structure behind modern single-cell analysis workflows.
Understanding its structure is essential for working with Scanpy, scverse, and large-scale transcriptomics datasets.
</p>

<ul>
  <li>Efficient data storage</li>
  <li>Flexible metadata handling</li>
  <li>Supports all major scRNA-seq analyses</li>
</ul>

---
