<h1 id="top">Getting Started with AnnData Package</h1>

<p>
A complete beginner-friendly guide to understanding and exploring AnnData objects in Python using Scanpy and scverse tools.
This tutorial demonstrates how single-cell RNA-seq data stored in <b>.h5ad</b> files can be loaded, inspected, and interpreted.
</p>
<h2>Table of Contents</h2>

<ul>
  <li><a href="#top">Getting Started with AnnData Package</a></li>

  <li><a href="#installation">Installation</a></li>

  <li><a href="#overview">Overview of AnnData Object</a></li>

  <li><a href="#first-look">First Look at AnnData Object</a></li>

  <li><a href="#active-matrix">Active Data Matrix (X)</a></li>

  <li><a href="#layers">Alternative Data in Layers</a></li>

  <li><a href="#obs-var">Cell and Gene Annotations (obs & var)</a></li>

  <li><a href="#indexing">Index Cells and Genes with obs_names and var_names</a></li>

  <li><a href="#subsetting-ann-data">Subsetting AnnData Objects</a></li>

  <li><a href="#subset-methods">Subsetting Methods</a></li>

  <li><a href="#basic-subsetting-example">Basic Subsetting Example</a></li>

  <li><a href="#boolean-subsetting">Boolean Subsetting Example</a></li>

  <li><a href="#multidimensional-annotations">Multidimensional Annotations (obsm & varm)</a></li>

  <li><a href="#obsp-example">Cell-Cell Distance Matrix (obsp)</a></li>

  <li><a href="#reordered-matrix">Reordered Distance Matrix</a></li>

  <li><a href="#uns-annotations">Unstructured Metadata (uns)</a></li>

  <li><a href="#views-copies">Views and Copies of AnnData Objects</a></li>

  <li><a href="#what-is-view">What is a View?</a></li>

  <li><a href="#creating-view">Creating a View</a></li>

  <li><a href="#view-shared-memory">Views Share Memory</a></li>

  <li><a href="#copy-vs-view">Copy vs View</a></li>

  <li><a href="#independent-objects">Independent Objects After Copy</a></li>

  <li><a href="#summary-views-copies">Summary (Views & Copies)</a></li>

  <li><a href="#summary">Final Summary</a></li>
</ul>


<h2 id="installation">Installation</h2>

<p>Install the required libraries before starting:</p>

<pre>
import anndata
import matplotlib.pyplot as plt
import numpy as np
import pooch
import scanpy as sc
</pre>

<h2 id="overview">Overview of AnnData Object</h2>

<p>
We use a preprocessed PBMC dataset hosted on figshare. It is downloaded using <b>pooch</b> and loaded into an AnnData object.
</p>

<pre>
datapath = pooch.retrieve(
    path=pooch.os_cache("scverse_tutorials"),
    url="https://exampledata.scverse.org/tutorials/scverse-getting-started-anndata-pbmc3k_processed.h5ad",
    known_hash="md5:b80deb0997f96b45d06f19c694e46243",
)

adata = anndata.read_h5ad(datapath)
</pre>

<h2 id="first-look">First Look at AnnData Object</h2>

<pre>
adata
</pre>

<pre>
AnnData object with n_obs × n_vars = 2638 × 11505
    obs: 'n_genes', 'percent_mito', 'n_counts', 'louvain_cell_types'
    var: 'gene_names', 'n_cells', 'gene_ids'
    uns: 'louvain', 'louvain_colors', 'pca'
    obsm: 'X_pca', 'X_tsne', 'X_umap'
    layers: 'raw'
    obsp: 'distances_all'
</pre>

<p>
This structure tells us:
</p>

<ul>
  <li>Number of cells and genes</li>
  <li>Cell-level metadata (obs)</li>
  <li>Gene-level metadata (var)</li>
  <li>Embeddings (obsm)</li>
  <li>Pairwise relationships (obsp)</li>
  <li>Additional metadata (uns)</li>
</ul>

<h2 id="active-matrix">Active Data Matrix (X)</h2>

<p>
The <b>adata.X</b> matrix contains the main expression data used for analysis. Rows represent cells and columns represent genes.
</p>

<pre>
adata.X
</pre>

<pre>
<Compressed Sparse Row sparse matrix of dtype 'float32'
 with 2076576 stored elements and shape (2638, 11505)>
</pre>

<h3>Sparse Matrix Properties</h3>

<pre>
print(adata.X.data)
print(adata.X.indices)
print(adata.X.nnz / np.prod(adata.X.shape))
</pre>

<p>
The data is stored in sparse format, meaning only non-zero values are saved to reduce memory usage.
</p>

<h2 id="layers">Alternative Data in Layers</h2>

<p>
Layers store alternative versions of the same dataset with identical shape as <b>adata.X</b>.
</p>

<pre>
adata.layers
</pre>

<pre>
Layers with keys: raw
</pre>

<h3>Raw Counts Layer</h3>

<pre>
adata.layers["raw"]
</pre>

<pre>
<Compressed Sparse Row sparse matrix of dtype 'int64'
 with 2076576 stored elements and shape (2638, 11505)>
</pre>

<h3>Create New Layer (CPM Normalization)</h3>

<pre>
adata.layers["counts_per_million"] = adata.layers["raw"].copy()

sc.pp.normalize_total(
    adata,
    target_sum=10**6,
    layer="counts_per_million"
)
</pre>

<pre>
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 2.5))
genes_of_interest = ["CD8A", "CD4", "KLRB1"]
sc.pl.matrixplot(
    adata,
    groupby="louvain_cell_types",
    var_names=genes_of_interest,
    layer="counts_per_million",  ## set which layer to plot
    ax=ax1,
    show=False,
)
ax1.set_title("CPM normalization")

sc.pl.matrixplot(
    adata,
    groupby="louvain_cell_types",
    var_names=genes_of_interest,
    layer="raw",  ## set which layer to plot
    ax=ax2,
    show=False,
)
ax2.set_title("raw counts")
plt.tight_layout()
</pre>

<img width="990" height="240" alt="normalization" src="https://github.com/user-attachments/assets/783b44a7-b783-4b9c-95c2-0b2e52c1d50c" />



<h2 id="obs-var">Cell and Gene Annotations (obs & var)</h2>

<p>
Metadata is stored in structured tables:
</p>

<ul>
  <li><b>obs:</b> cell-level annotations</li>
  <li><b>var:</b> gene-level annotations</li>
</ul>

<h3>Cell Metadata (obs)</h3>

<pre>
adata.obs
</pre>


<pre>
adata.obs.keys()
</pre>

<table border="1">
  <tr>
    <th>cell_barcode</th>
    <th>n_genes</th>
    <th>percent_mito</th>
    <th>n_counts</th>
    <th>louvain_cell_types</th>
  </tr>

  <tr>
    <td>AAACATACAACCAC-1</td>
    <td>781</td>
    <td>0.030178</td>
    <td>2419.0</td>
    <td>CD4 T cells</td>
  </tr>

  <tr>
    <td>AAACATTGAGCTAC-1</td>
    <td>1352</td>
    <td>0.037936</td>
    <td>4903.0</td>
    <td>B cells</td>
  </tr>

  <tr>
    <td>AAACATTGATCAGC-1</td>
    <td>1131</td>
    <td>0.008897</td>
    <td>3147.0</td>
    <td>CD4 T cells</td>
  </tr>

  <tr>
    <td>AAACCGTGCTTCCG-1</td>
    <td>960</td>
    <td>0.017431</td>
    <td>2639.0</td>
    <td>CD14+ Monocytes</td>
  </tr>

  <tr>
    <td>AAACCGTGTATGCG-1</td>
    <td>522</td>
    <td>0.012245</td>
    <td>980.0</td>
    <td>NK cells</td>
  </tr>

</table>
<p>As from this column we can see that we currently have four cell annotations:

i. the number of genes per cell

ii. fraction of counts from mitochondrial genes in each cell (a quality metric)

iii. total counts per cell

iv. cell type assignment of a louvain clustering

All of these annotations are one-dimensional & they contain a single number or string per cell. Also, each cell is uniquely identified by the index of the table. </p>

<pre>
Index(['n_genes', 'percent_mito', 'n_counts', 'louvain_cell_types'], dtype='object')
</pre>

<h3>Accessing a Column</h3>
<p> To access a specific column of the obs data frame, we can simply use its column name: </p>
<pre>
adata.obs["louvain_cell_types"]
</pre>
<h3>Cell Type Annotation (louvain_cell_types)</h3>

<table border="1">
  <tr>
    <th>cell_barcode</th>
    <th>louvain_cell_types</th>
  </tr>

  <tr>
    <td>AAACATACAACCAC-1</td>
    <td>CD4 T cells</td>
  </tr>

  <tr>
    <td>AAACATTGAGCTAC-1</td>
    <td>B cells</td>
  </tr>

  <tr>
    <td>AAACATTGATCAGC-1</td>
    <td>CD4 T cells</td>
  </tr>

  <tr>
    <td>AAACCGTGCTTCCG-1</td>
    <td>CD14+ Monocytes</td>
  </tr>

  <tr>
    <td>AAACCGTGTATGCG-1</td>
    <td>NK cells</td>
  </tr>

  <tr>
    <td>TTTCGAACTCTCAT-1</td>
    <td>CD14+ Monocytes</td>
  </tr>

  <tr>
    <td>TTTCTACTGAGGCA-1</td>
    <td>B cells</td>
  </tr>

  <tr>
    <td>TTTCTACTTCCTCG-1</td>
    <td>B cells</td>
  </tr>

  <tr>
    <td>TTTGCATGAGAGGC-1</td>
    <td>B cells</td>
  </tr>

  <tr>
    <td>TTTGCATGCCTCAC-1</td>
    <td>CD4 T cells</td>
  </tr>
</table>


<h3>Create New Annotation Column</h3>
<p> To add a new cell annotations, we simply add a new column to the obs data frame. For example, if we want to keep track of low-quality cells with a high mitochondrial gene count fraction, we can do the following to mark all cells above a certain threshold: </p>
<pre>
adata.obs["is_low_quality"] = adata.obs["percent_mito"] > 0.03
</pre>

<table border="1">
  <tr>
    <th>cell_barcode</th>
    <th>n_genes</th>
    <th>percent_mito</th>
    <th>n_counts</th>
    <th>louvain_cell_types</th>
    <th>is_low_quality</th>
  </tr>

  <tr>
    <td>AAACATACAACCAC-1</td>
    <td>781</td>
    <td>0.030178</td>
    <td>2419.0</td>
    <td>CD4 T cells</td>
    <td>True</td>
  </tr>

  <tr>
    <td>AAACATTGAGCTAC-1</td>
    <td>1352</td>
    <td>0.037936</td>
    <td>4903.0</td>
    <td>B cells</td>
    <td>True</td>
  </tr>

  <tr>
    <td>AAACATTGATCAGC-1</td>
    <td>1131</td>
    <td>0.008897</td>
    <td>3147.0</td>
    <td>CD4 T cells</td>
    <td>False</td>
  </tr>

  <tr>
    <td>AAACCGTGCTTCCG-1</td>
    <td>960</td>
    <td>0.017431</td>
    <td>2639.0</td>
    <td>CD14+ Monocytes</td>
    <td>False</td>
  </tr>

  <tr>
    <td>AAACCGTGTATGCG-1</td>
    <td>522</td>
    <td>0.012245</td>
    <td>980.0</td>
    <td>NK cells</td>
    <td>False</td>
  </tr>

  <tr>
    <td>TTTCGAACTCTCAT-1</td>
    <td>1155</td>
    <td>0.021104</td>
    <td>3459.0</td>
    <td>CD14+ Monocytes</td>
    <td>False</td>
  </tr>

  <tr>
    <td>TTTCTACTGAGGCA-1</td>
    <td>1227</td>
    <td>0.009294</td>
    <td>3443.0</td>
    <td>B cells</td>
    <td>False</td>
  </tr>

  <tr>
    <td>TTTCTACTTCCTCG-1</td>
    <td>622</td>
    <td>0.021971</td>
    <td>1684.0</td>
    <td>B cells</td>
    <td>False</td>
  </tr>

  <tr>
    <td>TTTGCATGAGAGGC-1</td>
    <td>454</td>
    <td>0.020548</td>
    <td>1022.0</td>
    <td>B cells</td>
    <td>False</td>
  </tr>

  <tr>
    <td>TTTGCATGCCTCAC-1</td>
    <td>724</td>
    <td>0.008065</td>
    <td>1984.0</td>
    <td>CD4 T cells</td>
    <td>False</td>
  </tr>

</table>

<h2 id="indexing">Index Cells and Genes with obs_names and var_names</h2>

<p>
A key feature of AnnData objects is shared indexing. Every cell and gene has a unique identifier that is consistently used across all data representations (X, obs, var, etc.).
This ensures that the same biological entity can always be referenced reliably.
</p>

<p>
Cell identifiers are stored in <code>obs_names</code>, and gene identifiers are stored in <code>var_names</code>.
</p>

<h3>Cell Index (obs_names)</h3>

<pre>
adata.obs_names
</pre>

<pre>
Index(['AAACATACAACCAC-1', 'AAACATTGAGCTAC-1', 'AAACATTGATCAGC-1',
       'AAACCGTGCTTCCG-1', 'AAACCGTGTATGCG-1', 'AAACGCACTGGTAC-1',
       'AAACGCTGACCAGT-1', 'AAACGCTGGTTCTT-1', 'AAACGCTGTAGCCA-1',
       'AAACGCTGTTTCTG-1',
       ...
       'TTTCAGTGTCACGA-1', 'TTTCAGTGTCTATC-1', 'TTTCAGTGTGCAGT-1',
       'TTTCCAGAGGTGAG-1', 'TTTCGAACACCTGA-1', 'TTTCGAACTCTCAT-1',
       'TTTCTACTGAGGCA-1', 'TTTCTACTTCCTCG-1', 'TTTGCATGAGAGGC-1',
       'TTTGCATGCCTCAC-1'],
      dtype='object', name='cell_barcode', length=2638)
</pre>

<h3>Gene Index (var_names)</h3>

<pre>
adata.var_names
</pre>

<pre>
Index(['LINC00115', 'NOC2L', 'KLHL17', 'PLEKHN1', 'HES4', 'ISG15', 'AGRN',
       'C1orf159', 'TNFRSF18', 'TNFRSF4',
       ...
       'MT-CO2', 'MT-ATP8', 'MT-ATP6', 'MT-CO3', 'MT-ND3', 'MT-ND4L', 'MT-ND4',
       'MT-ND5', 'MT-ND6', 'MT-CYB'],
      dtype='object', name='gene_names', length=11505)
</pre>

<hr>

<h2>Changing Gene Index</h2>

<p>
We can change how genes are indexed by assigning a different column from <code>adata.var</code> to <code>var_names</code>.
For example, we can use gene IDs instead of gene names.
</p>

<pre>
adata.var_names = adata.var["gene_ids"]
</pre>

<h3>Updated var Table</h3>

<pre>
adata.var
</pre>

<pre>
gene_names    n_cells    gene_ids
gene_ids
ENSG00000225880    LINC00115    18    ENSG00000225880
ENSG00000188976    NOC2L       258   ENSG00000188976
ENSG00000187961    KLHL17      9     ENSG00000187961
ENSG00000187583    PLEKHN1     7     ENSG00000187583
ENSG00000188290    HES4        145   ENSG00000188290
...
ENSG00000212907    MT-ND4L     398   ENSG00000212907
ENSG00000198886    MT-ND4      2588  ENSG00000198886
ENSG00000198786    MT-ND5      1399  ENSG00000198786
ENSG00000198695    MT-ND6      249   ENSG00000198695
ENSG00000198727    MT-CYB      2517  ENSG00000198727
</pre>

<p>
AnnData ensures that the length of the new index matches the number of genes. Otherwise, it raises an error.
</p>

<hr>

<h2>Restoring Gene Names</h2>

<p>
For better interpretability, we switch back to gene names as indices.
</p>

<pre>
adata.var_names = adata.var["gene_names"]
</pre>

<h3>Final var Table</h3>

<pre>
adata.var
</pre>

<pre>
gene_names    n_cells    gene_ids
gene_names
LINC00115    LINC00115    18    ENSG00000225880
NOC2L        NOC2L        258   ENSG00000188976
KLHL17       KLHL17       9     ENSG00000187961
PLEKHN1      PLEKHN1      7     ENSG00000187583
HES4         HES4         145   ENSG00000188290
...
MT-ND4L      MT-ND4L      398   ENSG00000212907
MT-ND4       MT-ND4       2588  ENSG00000198886
MT-ND5       MT-ND5       1399  ENSG00000198786
MT-ND6       MT-ND6       249   ENSG00000198695
MT-CYB       MT-CYB       2517  ENSG00000198727
</pre>
<h2 id="subsetting-ann-data">Subsetting AnnData Objects</h2>

<p>
AnnData objects support flexible subsetting of cells (observations) and genes (variables) using consistent indexing rules.
This allows efficient filtering and exploration of single-cell datasets without breaking the data structure.
</p>

<h3 id="subset-methods">Subsetting Methods</h3>

<ul>
  <li><b>Name-based indexing:</b> Using <code>obs_names</code> or <code>var_names</code> (e.g., gene names or cell barcodes)</li>
  <li><b>Numerical indexing:</b> Using integer slices such as <code>:10</code>, <code>0:10</code>, or lists like <code>[0, 2, 4]</code></li>
  <li><b>Boolean indexing:</b> Using conditions like <code>adata.obs["is_low_quality"]</code> or cluster-based filtering</li>
</ul>

<p>
All subsetting operations automatically update:
</p>

<ul>
  <li><code>adata.X</code> (active matrix)</li>
  <li>All <code>layers</code></li>
  <li><code>obs</code> and <code>var</code> annotations</li>
</ul>

<h3 id="basic-subsetting-example">Basic Subsetting Example</h3>

<pre><code class="language-python">
adata_small = adata[:5, ["LYZ", "FOS", "MALAT1"]]
adata_small.shape
</code></pre>

<pre>
(5, 3)
</pre>

<h4>Active Data Matrix</h4>

<pre><code class="language-python">
adata_small.X.toarray()
</code></pre>

<pre>
array([[0.6496621 , 2.4036813 , 3.8249342 ],
       [0.8553989 , 0.64266497, 4.174533  ],
       [0.878057  , 2.8014445 , 4.7918878 ],
       [3.0494576 , 1.6464083 , 2.3237872 ],
       [0.        , 0.        , 3.923651  ]], dtype=float32)
</pre>

<h4>Raw Layer After Subsetting</h4>

<pre><code class="language-python">
adata_small.layers["raw"].toarray()
</code></pre>

<pre>
array([[  1,  11,  49],
       [  3,   2, 142],
       [  2,  22, 170],
       [ 24,   5,  11],
       [  0,   0,  22]])
</pre>

<h4>Cell Metadata (obs)</h4>

<pre><code class="language-python">
adata_small.obs
</code></pre>

<pre>
n_genes  percent_mito  n_counts  louvain_cell_types   is_low_quality
AAACATACAACCAC-1   781   0.030178   2419.0   CD4 T cells        True
AAACATTGAGCTAC-1  1352   0.037936   4903.0   B cells            True
AAACATTGATCAGC-1  1131   0.008897   3147.0   CD4 T cells        False
AAACCGTGCTTCCG-1   960   0.017431   2639.0   CD14+ Monocytes    False
AAACCGTGTATGCG-1   522   0.012245    980.0   NK cells           False
</pre>

<h4>Gene Metadata (var)</h4>

<pre><code class="language-python">
adata_small.var
</code></pre>

<pre>
gene_names   n_cells   gene_ids
LYZ          LYZ       1631      ENSG00000090382
FOS          FOS       2473      ENSG00000170345
MALAT1       MALAT1    2699      ENSG00000251562
</pre>

<h3 id="boolean-subsetting">Boolean Subsetting Example</h3>

<pre><code class="language-python">
adata_high_quality = adata[~adata.obs["is_low_quality"], :]
adata_high_quality.obs
</code></pre>

<pre>
2257 rows × 5 columns
</pre>

<p>
Boolean indexing is commonly used to remove low-quality cells or filter specific biological populations.
</p>

<hr>

<h2 id="multidimensional-annotations">Multidimensional Annotations for Cells and Genes (obsm & varm)</h2>

<p>
Sometimes, simple one-dimensional annotations in <code>obs</code> and <code>var</code> are not enough to store complex results of analysis.
For example, dimensionality reduction methods like PCA, t-SNE, and UMAP generate multi-dimensional representations of cells.
</p>

<p>
To store such structured data, AnnData provides:
</p>

<ul>
  <li><code>obsm</code> → multidimensional annotations for cells</li>
  <li><code>varm</code> → multidimensional annotations for genes</li>
</ul>

<h3>Available Keys in uns</h3>

<pre>
adata.uns.keys()
</pre>

<pre>
dict_keys(['louvain', 'louvain_colors', 'pca'])
</pre>

<h3>Cell-Level Multidimensional Data (obsm)</h3>

<pre>
adata.obsm
</pre>

<pre>
AxisArrays with keys: X_pca, X_tsne, X_umap
</pre>

<p>
These represent different embeddings of cells stored as matrices inside a dictionary-like structure.
</p>

<h3>Shape of Each Embedding</h3>

<pre>
for key in adata.obsm:
    print(key, adata.obsm[key].shape)
</pre>

<pre>
X_pca   (2638, 50)
X_tsne  (2638, 2)
X_umap  (2638, 2)
</pre>

<p>
Each row corresponds to a cell, while columns represent reduced dimensions.
The number of cells always matches <code>adata.X</code>, but dimensionality depends on the method.
</p>

<hr>

<h3>Visualization of Embeddings</h3>

<p>
We can visualize PCA, t-SNE, and UMAP embeddings and highlight B cells for interpretation.
</p>

<h4>PCA Plot</h4>

<pre>
plt.scatter(
    x=adata.obsm["X_pca"][:, 0],
    y=adata.obsm["X_pca"][:, 1],
    c=adata.obs["louvain_cell_types"] == "B cells",
    s=3
)
plt.title("PCA")
</pre>

<h4>t-SNE Plot</h4>

<pre>
plt.scatter(
    x=adata.obsm["X_tsne"][:, 0],
    y=adata.obsm["X_tsne"][:, 1],
    c=adata.obs["louvain_cell_types"] == "B cells",
    s=3
)
plt.title("t-SNE")
</pre>

<h4>UMAP Plot</h4>

<pre>
plt.scatter(
    x=adata.obsm["X_umap"][:, 0],
    y=adata.obsm["X_umap"][:, 1],
    c=adata.obs["louvain_cell_types"] == "B cells",
    s=3
)
plt.title("UMAP")
</pre>
<img width="945" height="350" alt="embedded visualization" src="https://github.com/user-attachments/assets/ef99e73a-f508-4fa8-9a73-e75eabc09974" />



<h3>Key Insight</h3>

<p>
<code>obsm</code> allows storage of multi-dimensional cell representations, while <code>varm</code> performs the same role for genes.
These structures are essential for storing PCA embeddings, UMAP/t-SNE coordinates, and gene loadings in a consistent way.
</p>

<p>
This design keeps AnnData flexible for modern single-cell analysis workflows.
</p>

<h2 id="pairwise-annotations">Pairwise Annotations (obsp & varp)</h2>

<p>
Some analyses require relationships between pairs of cells or genes, such as:
</p>

<ul>
  <li>Cell-cell distance matrices</li>
  <li>K-nearest neighbor graphs</li>
  <li>Gene-gene interaction networks</li>
</ul>

<p>
These are stored in:
</p>

<ul>
  <li><code>obsp</code> → cell-cell matrices</li>
  <li><code>varp</code> → gene-gene matrices</li>
</ul>

<h3 id="obsp-example">Cell-Cell Distance Matrix</h3>

<pre><code class="language-python">
adata.obsp["distances_all"]
</code></pre>

<pre>
(n_cells × n_cells matrix)
</pre>

<h4>Distance Matrix Visualization</h4>

<pre><code class="language-python">
plt.imshow(adata.obsp["distances_all"])
plt.colorbar()
plt.show()
</code></pre>
<p align="center">
<img width="524" height="418" alt="distance matrix visualization" src="https://github.com/user-attachments/assets/71db8b14-1091-43db-b7ef-e99f9772137d" />
</p>


<h3 id="reordered-matrix">Reordered Distance Matrix</h3>

<pre><code class="language-python">
reorder = np.argsort(adata.obs["louvain_cell_types"])
plt.imshow(adata[reorder, :].obsp["distances_all"])
plt.colorbar()
plt.show()
</code></pre>
<p align="center">
  <img width="524" height="418" alt="Reorder Distance Matrix" src="https://github.com/user-attachments/assets/2f2a475f-bfbc-4977-8dde-d21e2fc5d690" />
</p>

<p>
Reordering cells by cluster reveals clear block structures in the distance matrix.
</p>

<hr>

<h2 id="uns-annotations">Unstructured Metadata (uns)</h2>

<p>
The <code>uns</code> field stores dataset-level metadata that does not belong to specific cells or genes.
It is commonly used for:
</p>

<ul>
  <li>Clustering parameters</li>
  <li>Color palettes</li>
  <li>PCA variance information</li>
  <li>Method outputs and configurations</li>
</ul>

<h3 id="uns-structure">Structure of uns</h3>

<pre><code class="language-python">
adata.uns.keys()
</code></pre>

<pre>
dict_keys(['louvain', 'louvain_colors', 'pca'])
</pre>

<h3 id="uns-example">Example Contents</h3>

<pre><code class="language-python">
adata.uns["louvain"]
</code></pre>

<pre>
{'params': {'random_state': array([0]), 'resolution': array([1])}}
</pre>

<pre><code class="language-python">
adata.uns["pca"]["variance"]
</code></pre>

<pre>
Explained variance values for principal components
</pre>

<p>
The <code>uns</code> field ensures that dataset-wide metadata is preserved throughout analysis workflows.
</p>

<h2 id="views-copies">Views and Copies of AnnData Objects</h2>

<p>
In AnnData, subsetting behavior is designed to be memory efficient. When you subset an AnnData object, you typically get a <b>view</b> instead of a full copy. This avoids duplicating large expression matrices in memory, which is critical for large-scale single-cell datasets.
</p>

<p>
However, in some cases you may want an independent object. For that, you explicitly create a <b>copy</b>. Understanding the difference between views and copies is essential for safe and efficient data manipulation.
</p>

<h3 id="what-is-view">What is a View?</h3>

<p>
A view is a lightweight reference to the original AnnData object. It does not store its own data; instead, it points back to the parent object. Any changes in the parent object are reflected in the view.
</p>

<pre><code class="language-python">
adata.X[:5, 5:10].toarray()
</code></pre>

<pre>
array([[0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.64266497, 0.        ],
       [0.532456  , 0.        , 0.        , 0.        , 0.        ],
       [2.1446393 , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ]],
      dtype=float32)
</pre>

<h3 id="creating-view">Creating a View</h3>

<p>
When you subset AnnData directly, you get a view instead of a copy:
</p>

<pre><code class="language-python">
adata_view = adata[:5, 5:10]
adata_view
</code></pre>

<pre>
View of AnnData object with n_obs × n_vars = 5 × 5
    obs: 'n_genes', 'percent_mito', 'n_counts', 'louvain_cell_types', 'is_low_quality'
    var: 'gene_names', 'n_cells', 'gene_ids'
    uns: 'louvain', 'louvain_colors', 'pca'
    obsm: 'X_pca', 'X_tsne', 'X_umap'
    layers: 'raw', 'counts_per_million'
    obsp: 'distances_all'
</pre>

<h3 id="view-shared-memory">Views Share Memory</h3>

<p>
Since views are linked to the parent object, modifying the parent will also affect the view.
</p>

<pre><code class="language-python">
adata.X[0, 7] = 99
adata_view.X.toarray()
</code></pre>

<pre>
array([[ 0.        ,  0.        , 99.        ,  0.        ,  0.        ],
       [ 0.        ,  0.        ,  0.        ,  0.64266497,  0.        ],
       [ 0.532456  ,  0.        ,  0.        ,  0.        ,  0.        ],
       [ 2.1446393 ,  0.        ,  0.        ,  0.        ,  0.        ],
       [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ]],
      dtype=float32)
</pre>

<p>
This confirms that the view directly reflects changes in the parent AnnData object.
</p>

<div class="note">
<b>Note:</b> Propagation of changes from parent to view works for numpy-based fields like X, layers, obsm, varm, obsp, and varp.
However, it does NOT work for pandas-based fields such as obs and var.
</div>

<h3 id="copy-vs-view">Turning a View into a Copy</h3>

<p>
A view can be converted into an independent AnnData object in two ways:
</p>

<ul>
  <li>By calling <code>.copy()</code></li>
  <li>By modifying the view (which triggers an implicit copy)</li>
</ul>

<h4>Using .copy()</h4>

<pre><code class="language-python">
adata_view_copy = adata_view.copy()
adata_view_copy
</code></pre>

<pre>
AnnData object with n_obs × n_vars = 5 × 5
    obs: 'n_genes', 'percent_mito', 'n_counts', 'louvain_cell_types', 'is_low_quality'
    var: 'gene_names', 'n_cells', 'gene_ids'
    uns: 'louvain', 'louvain_colors', 'pca'
    obsm: 'X_pca', 'X_tsne', 'X_umap'
    layers: 'raw', 'counts_per_million'
    obsp: 'distances_all'
</pre>

<h4>Implicit Conversion to Copy</h4>

<pre><code class="language-python">
adata_view.obs["new_column"] = "Test"
adata_view
</code></pre>

<pre>
AnnData object with n_obs × n_vars = 5 × 5
    obs: 'n_genes', 'percent_mito', 'n_counts', 'louvain_cell_types', 'is_low_quality', 'new_column'
    var: 'gene_names', 'n_cells', 'gene_ids'
    uns: 'louvain', 'louvain_colors', 'pca'
    obsm: 'X_pca', 'X_tsne', 'X_umap'
    layers: 'raw', 'counts_per_million'
    obsp: 'distances_all'
</pre>

<p>
At this point, the view is converted into a full AnnData object internally.
</p>

<h3 id="independent-objects">Independent Objects After Copy</h3>

<p>
Once a view is converted into a copy, it becomes independent from the parent object.
Changes in the original data will no longer affect it.
</p>

<pre><code class="language-python">
adata.X[0, 8] = 98
adata.X[:5, 5:10].toarray()
</code></pre>

<pre>
array([[ 0.        ,  0.        , 99.        , 98.        ,  0.        ],
       [ 0.        ,  0.        ,  0.        ,  0.64266497,  0.        ],
       [ 0.532456  ,  0.        ,  0.        ,  0.        ,  0.        ],
       [ 2.1446393 ,  0.        ,  0.        ,  0.        ,  0.        ],
       [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ]],
      dtype=float32)
</pre>

<pre><code class="language-python">
adata_view.X.toarray()
adata_view_copy.X.toarray()
</code></pre>

<p>
Both the converted view and its copy remain unchanged by further modifications in the parent object.
</p>

<h3 id="summary-views-copies">Summary</h3>

<ul>
  <li>Views are memory-efficient and reference the original AnnData object.</li>
  <li>Copies are independent objects with their own stored data.</li>
  <li>Modifying a view may trigger an automatic conversion into a full AnnData object.</li>
  <li>Use views for analysis workflows; use copies when you need independent manipulation.</li>
</ul>
<h2 id="summary">Summary</h2>
<p>
AnnData combines all single-cell data components into one structured object.
</p>

<ul>
  <li>X → Expression matrix</li>
  <li>obs → Cell metadata</li>
  <li>var → Gene metadata</li>
  <li>obsm → Embeddings</li>
  <li>obsp → Pairwise relationships</li>
  <li>uns → Global metadata</li>
</ul>
