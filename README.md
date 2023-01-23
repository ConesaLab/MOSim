# MOSim
MOSim is an R package for the simulation of multi-omic experiments that mimic regulatory mechanisms within the cell. 
Gene expression (RNA-seq count data) is the central data type simulated by MOSim, while the rest of available omic data types 
provide gene regulation information and include ATAC-seq (DNase-seq), ChIP-seq, small RNA-seq and Methyl-seq. In addition to these omics, 
regulation by transcription factors (TFs) can also be modeled.

### Installation

MOSim is a Bioconductor R package, and we strongly recommend that it is installed from the Bioconductor repository. 
To install MOSim, open the R console and run:

  ```
  if (!requireNamespace("BiocManager", quietly = TRUE)) 
      install.packages("BiocManager")

  BiocManager::install("MOSim")
  ```

The developer version can be installed from GitHub using the devtools R package:
	
  ```
  install.packages("devtools")
  devtools::install_github("ConesaLab/MOSim")
  ```

### Documentation

Vignettes and documentation can be accessed from [MOSim's Bioconductor site](http://bioconductor.org/packages/release/bioc/html/MOSim.html), 
or by running the following line in the R console:

	browseVignettes("MOSim")


### Citation

If you used MOSim for your research, please cite:

- Martínez C, Tarazona S (2020). MOSim: Multi-Omics Simulation (MOSim). R package version 1.4.0, https://github.com/ConesaLab/MOSim.
