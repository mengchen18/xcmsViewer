# xcmsViewer

The xcmsViewer is an R packages used for analysing and visualising the 
untargeted metabolomics data. It emphasizes the following features:

* **Annotation/metabolite identification**: The feature annotation (or identification of metabolite)
are performed on both MS1 and MS2 levels. On the MS1 level, the mass of metabolites (by adding or 
removing certain adducts) are compared with the monoisotopic mass of known metabolites to annotate features. 
On the MS2 level, the MS2 spectrum is compared with the spectrum of know metabolites. The information 
could be downloaded from public databases, such as HMDB or MS-DIAL, or users' own databases.

* **Expression matrix-centric**: The feature table given by xcms package is a matrix, where the rows 
are features/metabolites, the columns are samples measured and the values in the matrix are the 
intensities of features/metabolites. Such a matrix is a starting point for downstream analysis, 
such as differential expression, feature selection, etc. In xcmsViewer, the "expression matrix" is 
the center of all other information (chromatogram, annotation using mass, annotation using MS2 
spectra, etc.

* **Interactive visualization**: The process of converting raw files to the expression matrix is rather 
well established in transcriptomics and proteomics. However, due to the immaturity of metabolomics 
data processing, we often need to check the extract chromatograms behind intensity and how well a
MS2 spectrum matched to a reference spectrum. The xcmsViewer offers interactive visualization of 
metabolite intensities, chromatograms and MS2 spectra.

* **Flexibility**: The xcmsViewer is seamlessly integrated into the R environment, therefore, it inherits 
the flexibility of statistical analysis given by R. 

* **Deliverable**: The processed data, together with the visualization using Shiny app, could be 
delievered to the collaborator as standalone software or using shiny-server. Your collaborator don't 
need to install any program to explore the results.


# Installation and examples
## xcmsViewerApp
Installation
```{r}
devtools::install_github("mengchen18/xcmsViewer")
```
This command will open a shiny app. User can select the dataset to load (there is only one as an example) using the selection box on the topright corner.
