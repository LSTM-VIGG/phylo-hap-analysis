# IQ-TREE: A Beginner's Guide

IQ-TREE is a powerful, fast, and user-friendly software for phylogenetic analysis using maximum likelihood methods. This guide provides an introduction to IQ-TREE concepts and usage, specifically in the context of haplotype phylogenetic analysis.

## What is IQ-TREE?

IQ-TREE is a state-of-the-art software package for reconstructing phylogenetic trees from molecular sequence data. It's designed to be both fast and accurate, making it ideal for large-scale genomic analyses.

### Key Features

- **Fast computation**: Uses advanced algorithms and multi-threading
- **Automatic model selection**: ModelFinder Plus automatically chooses the best evolutionary model
- **Bootstrap support**: Provides statistical confidence measures for tree branches
- **Large dataset handling**: Efficiently processes thousands of sequences
- **User-friendly**: Simple command-line interface with sensible defaults

## Basic Concepts

### What is Phylogenetic Reconstruction?

Phylogenetic reconstruction involves inferring evolutionary relationships between organisms (or in our case, haplotypes) based on their molecular sequences. The result is a tree-like diagram showing how sequences are related through evolutionary time.

### Maximum Likelihood Method

IQ-TREE uses the **maximum likelihood (ML)** approach, which:
- Finds the tree topology and branch lengths that best explain the observed sequence data
- Uses statistical models of molecular evolution
- Provides a principled framework for comparing different evolutionary hypotheses

### Key Terms

- **Tree topology**: The branching pattern showing relationships
- **Branch lengths**: Represent evolutionary distance (number of substitutions)
- **Bootstrap values**: Statistical confidence measures (0-100%)
- **Outgroup**: A reference sequence used to root the tree
- **Model of evolution**: Mathematical description of how sequences change over time

## Understanding Evolutionary Models

### Why Models Matter

Different genes and genomic regions evolve at different rates and with different patterns. IQ-TREE uses mathematical models to account for:
- Different substitution rates between nucleotides
- Rate variation across sites
- Base composition bias

### ModelFinder Plus (MFP)

Our pipeline uses `MFP` (ModelFinder Plus), which:
- Tests hundreds of evolutionary models
- Selects the best-fitting model using statistical criteria
- Accounts for both model complexity and fit to data
- Saves time by not requiring manual model selection

### Common Models You Might See

- **GTR+F+I+G**: General time-reversible model with empirical base frequencies, invariant sites, and gamma-distributed rate variation
- **HKY+F+G**: Hasegawa-Kishino-Yano model with empirical base frequencies and gamma rates
- **JC**: Jukes-Cantor model (simplest, assumes equal rates and frequencies)

## Understanding IQ-TREE Output

### File Extensions and Their Meanings

When IQ-TREE runs on your FASTA file (e.g., `vgsc_focal.fasta`), it creates several output files:

#### `.treefile` - The Final Tree
- Contains the best ML tree in Newick format
- This is what you'll use for visualization and downstream analysis
- Branch lengths represent expected substitutions per site

#### `.iqtree` - Main Results File
Contains detailed analysis results:
```
Model of evolution: GTR+F+I+G4
Log-likelihood: -12345.67
Tree length: 0.12345
```

#### `.log` - Analysis Log
- Complete record of the analysis
- Shows model testing results, optimization progress
- Useful for troubleshooting

#### `.contree` - Consensus Tree
- Bootstrap consensus tree (if bootstrapping was performed)
- Shows support values as percentages

### Reading Tree Files

Trees are stored in **Newick format**:
```
((Sample1:0.001,Sample2:0.002):0.005,Sample3:0.008);
```

- Parentheses group related sequences
- Colons precede branch lengths
- Semicolon ends the tree
- Numbers after colons are evolutionary distances

### Interpreting Bootstrap Values

Bootstrap values indicate statistical support for branches:
- **95-100%**: Very strong support
- **85-94%**: Strong support  
- **75-84%**: Moderate support
- **<75%**: Weak support

## Command Line Usage

### Basic IQ-TREE Command
```bash
iqtree -s alignment.fasta -m MFP -bb 1000
```

Where:
- `-s`: Input alignment file
- `-m MFP`: Use ModelFinder Plus for model selection
- `-bb 1000`: Perform 1000 bootstrap replicates

### Our Pipeline's Command
The `run_iqtree.sh` script uses:
```bash
iqtree2 -s input.fasta -m MFP -bb 1000 -T AUTO
```

Additional parameters:
- `-T AUTO`: Automatically determine optimal number of CPU threads
- `-o outgroup`: Specify outgroup for rooting (optional)
- `--prefix name`: Set custom prefix for output files

### Common Options

```bash
# Basic analysis with model selection
iqtree -s data.fasta -m MFP

# With bootstrap support
iqtree -s data.fasta -m MFP -bb 1000

# Specify number of threads
iqtree -s data.fasta -m MFP -T 4

# Set outgroup
iqtree -s data.fasta -m MFP -o OutgroupName

# Custom output prefix
iqtree -s data.fasta -m MFP --prefix my_analysis
```

## Tips for Better Results

### Input Data Quality
- **Remove gaps**: IQ-TREE handles gaps, but excessive gaps can affect accuracy
- **Check alignment**: Ensure sequences are properly aligned
- **Adequate sampling**: More sequences generally improve accuracy
- **Informative sites**: Ensure sufficient variable positions

### Performance Optimization
- **Use multiple cores**: `-T AUTO` or `-T 8` for 8 cores
- **Bootstrap wisely**: 1000 replicates is usually sufficient
- **Consider memory**: Large datasets may need memory management options

### Model Selection
- **Trust ModelFinder**: Let `MFP` choose the best model
- **Don't over-parameterize**: Complex models aren't always better
- **Check model adequacy**: Review the `.iqtree` file for warnings

## Interpreting Results in Context

### For Haplotype Analysis

When analyzing haplotypes (as in this pipeline):
- **Expect short branch lengths**: Haplotypes from the same species are closely related
- **Look for clustering**: Samples from same populations/species should group together
- **Consider recombination**: Phylogenetic signal may be limited in recombining regions
- **Use multiple regions**: Compare focal, upstream, and downstream analyses

### Biological Interpretation

Strong bootstrap support suggests:
- **Reliable relationships**: High confidence in branching pattern
- **Consistent signal**: Multiple sites support the same topology
- **Evolutionary significance**: Likely reflects true evolutionary history

Weak bootstrap support may indicate:
- **Rapid evolution**: Not enough time for mutations to accumulate
- **Recombination**: Conflicting evolutionary histories
- **Insufficient data**: Need more sequence information

## Further Reading

### Essential Papers
- **Nguyen et al. (2015)**: Original IQ-TREE paper in Mol Biol Evol
- **Kalyaanamoorthy et al. (2017)**: ModelFinder paper in Nat Methods
- **Hoang et al. (2018)**: UFBoot2 bootstrap method in Mol Biol Evol

### Online Resources
- [IQ-TREE Website](http://www.iqtree.org/): Official documentation and tutorials
- [IQ-TREE Manual](http://www.iqtree.org/doc/): Comprehensive usage guide
- [IQ-TREE Tutorial](http://www.iqtree.org/workshop/): Workshop materials