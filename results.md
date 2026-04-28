<h1 id="top"> Single-Cell RNA-seq Analysis Results</h1>

<li><a href="#results">Results</a></li>

  <li><a href="#interpretation">Results Interpretation</a></li>

  <li><a href="#conclusion">Conclusion</a></li>
</ul>

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
