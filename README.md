# Phylo Hap Analysis 

Code for phylogenetic analysis of haplotypes from malariagen_data. It downloads haplotype data from the cloud, generates multiple sequence alignments in FASTA format, constructs phylogenetic trees using IQ-TREE, and provides interactive visualizations using either ape (R) or baltic_bokeh (Python).

## Overview

The notebooks/scripts perform the following steps:
1. **Data extraction**: Downloads haplotypes and metadata from MalariaGEN using the `malariagen_data` Python package
2. **FASTA generation**: Converts haplotype data to multiple sequence alignment format
3. **Phylogenetic reconstruction**: Uses IQ-TREE to infer maximum likelihood phylogenetic trees
4. **Visualization**: Creates plots using R (ape package) or Python (baltic_bokeh)

## Prerequisites

### System Requirements
- Python 3.8+
- R 4.0+
- IQ-TREE v2.0+ (phylogenetic inference)

## Installation

1. Clone this repository:
```bash
git clone https://github.com/LSTM-VIGG/phylo-hap-analysis
cd phylo-hap-analysis
```

2. Install Python dependencies:
```bash
pip install malariagen_data jupyter
```

3. Install baltic_bokeh from GitHub:
```bash
git clone https://github.com/sanjaynagi/baltic-bokeh.git
cd baltic-bokeh
pip install -e .
cd ..
```

4. Install R packages (in R console):
```r
install.packages(c("ape", "phytools", "stringr", "data.table", "dplyr"))
```

5. Install IQ-TREE following [official instructions](http://www.iqtree.org/doc/Quickstart)

## Components

### 1. Data Extraction and FASTA Generation
**File**: `notebooks/write_phylo_hap_fastas.ipynb`

This Jupyter notebook:
- Connects to the MalariaGEN Ag3 API
- Downloads haplotype data for specified genomic regions
- Converts haplotypes to FASTA format for phylogenetic analysis
- Generates metadata files with sample information
- Creates focal, upstream, and downstream region analyses

**Key parameters to configure**:
- `region`: Genomic coordinates of interest
- `locus_name`: Name for output files
- `flanking`: Flanking region size (bp)
- `sample_sets`: MalariaGEN cohorts to include

**Outputs**:
- `fastas/{locus_name}_{region}.fasta`: Multiple sequence alignments
- `fastas/{locus_name}_{region}.metadata.tsv`: Sample metadata

### 2. Phylogenetic Tree Inference
**File**: `scripts/run_iqtree.sh`

Bash script that runs IQ-TREE on all FASTA files:
- Automatically processes all `.fasta` files in the `fastas/` directory
- Uses ModelFinder Plus (MFP) for automatic model selection
- Performs bootstrap analysis (default: 1000 replicates)
- Supports automatic threading

**Configuration**:
```bash
FASTA_DIR="fastas"                    # Input directory
OUTPUT_DIR="iqtree_results"           # Output directory
OUTGROUP="VBS00259-4651STDY7017186_1" # Outgroup sample ID
BOOTSTRAPS=1000                       # Bootstrap replicates
THREADS="AUTO"                        # CPU threads
```

**Usage**:
```bash
chmod +x scripts/run_iqtree.sh
./scripts/run_iqtree.sh
```

### 3. Tree Visualization

#### Option A: R Visualization (APE)
**Files**: `scripts/xavi-ape-plotting.R`, `scripts/ape-plotting.R`

R scripts using the `ape` package for phylogenetic plotting:
- Reads IQ-TREE output files (`.treefile`)
- Applies midpoint rooting
- Color-codes tips by species, karyotype, or other metadata
- Generates PDF outputs with publication-ready figures

**Key features**:
- Automatic tip coloring by taxonomic groups
- Branch length constraints for visualization
- Multiple output formats (standard and large format)

**Usage**:
```r
# Edit file paths in the script, then run:
Rscript scripts/ape-plotting.R
```

#### Option B: Interactive Python Visualization (baltic_bokeh)
**File**: `scripts/bokeh-phy.py`

Interactive web-based phylogenetic visualization:
- Bokeh-based interactive plots
- Hover tooltips with metadata
- Customizable color schemes
- WebGL rendering for large trees

**Usage** (after conversion to command-line script):
```bash
python scripts/bokeh-phy.py --tree results/tree.newick --metadata results/metadata.tsv --output plot.html
```

## Example

1. **Configure parameters** in `notebooks/write_phylo_hap_fastas.ipynb`
2. **Run the notebook** to generate FASTA files
3. **Run phylogenetic inference**:
   ```bash
   chmod +x scripts/run_iqtree.sh
   ./scripts/run_iqtree.sh
   ```
4. **Visualize results**:
   - R: `Rscript scripts/ape-plotting.R`
   - Python: `python scripts/bokeh-phy.py [arguments]`

## Output Structure

```
phylo-hap-analysis/
├── fastas/                          # FASTA alignments and metadata
│   ├── vgsc_focal.fasta
│   ├── vgsc_focal.metadata.tsv
│   ├── vgsc_upstream.fasta
│   └── ...
├── iqtree_results/                  # IQ-TREE outputs
│   ├── vgsc_focal.fasta.treefile    # Final phylogenetic tree
│   ├── vgsc_focal.fasta.log         # Analysis log
│   ├── vgsc_focal.fasta.iqtree      # Detailed results
│   └── ...
├── notebooks/                       # Jupyter notebooks
├── scripts/                         # Analysis scripts
└── results_cache/                   # MalariaGEN data cache
```

## Support

- **IQ-TREE**: Consult [IQ-TREE manual](http://www.iqtree.org/doc/)
- [IQTREE.md](IQTREE.md) - Beginner's guide to IQ-TREE
- [MalariaGEN documentation](https://malariagen.github.io/malariagen-data-python/)
- [APE package tutorial](http://ape-package.ird.fr/ape_tutorial.pdf)
