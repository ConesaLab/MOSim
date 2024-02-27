# MOSim
MOSim is an R package for the simulation of multi-omic bulk and single cell experiments that mimic regulatory mechanisms within the cell. 
Gene expression (RNA-seq count data) is the central data type simulated by MOSim, while the rest of available omic data types 
provide gene regulation information. For bulk simulation, regulators include ATAC-seq (DNase-seq), ChIP-seq, miRNA-seq and Methyl-seq. In addition to these omics, 
regulation by transcription factors (TFs) can also be modeled. While for single-cell simulation, the regulator included is scATAC-seq.

### Installation

MOSim is a Bioconductor R package, and we strongly recommend that it is installed from the Bioconductor repository. 
To install MOSim, open the R console and run:

  ```
  if (!requireNamespace("BiocManager", quietly = TRUE)) 
      install.packages("BiocManager")

  BiocManager::install("MOSim")
  ```

The developer version (which now includes the scMOSim functionalities) can be installed from GitHub using the devtools R package:
	
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

- Monzó C, Martínez C, Arzalluz-Luque A,  Conesa A, Tarazona S (2024). MOSim: bulk and single-cell multi-layer regulatory network simulator. BioRxiv. DOI: 10.1101/421834
  
scMOSim strongly relies on functionality from SPARSim. If you use the scMOSim module to simulate multi-omics single cell data, please also cite:  
- Baruzzo G, Patuzzi I, Di Camillo B (2020). SPARSim single cell: a count data simulator for scRNA-seq data. Bioinformatics, Volume 36, Issue 5, Pages 1468-1475. DOI: 10.1093/bioinformatics/btz752
