# Oslo 2023: vsearch-swarm pipeline

[![Open In
Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/frederic-mahe/BIO9905MERG1_vsearch_swarm_pipeline/blob/main/pipeline.ipynb)

A [vsearch](https://github.com/torognes/vsearch),
[cutadapt](https://cutadapt.readthedocs.io/en/stable/) and
[swarm](https://github.com/torognes/swarm)-based pipeline for
metabarcoding data.

## edit and export

Edit with [Emacs](https://www.gnu.org/software/emacs/) in
[markdown](https://jupyter-notebook.readthedocs.io/en/stable/examples/Notebook/Working%20With%20Markdown%20Cells.html)
format, export to a [jupyter](https://jupyter.org/) notebook with
[pandoc](https://pandoc.org):

``` bash
pandoc -o pipeline.ipynb pipeline.md
```

## Tasks

- [x] check for typos and spelling,
- [x] export to `.ipynb`,
- [x] run on Google colab
