# COMBINe: Cell detectiOn in Mouse BraIN
Designed for automated detection of neurons and astrocytes in cleared MADM mouse brains.  
Datasets used in this project were acuiqred using a custom-built light sheet fluorescence microscope. 

The workflow contains the following steps:
1. Stitching* - TeraStitcher  (advanced mode recommended)
2. Channel alignment* (optional but recommended for light sheet datasets) - NuMorph (MATLAB)
3. Cell detection - fizyr/keras-retinanet (Python)  
4. Cell stitching (Python)  
5. Registration - ClearMap (Python)  
6. Analysis (MATLAB)  

*The order of these two steps doesn't matter.  

## Running COMBINe  

![The flow chart of running the pipeline](/image/flowchart.png "Running COMBINe")  
The pipeline is designed for image dataset/3D volume stored as [a two-level hierarchy of folders](https://github.com/abria/TeraStitcher/wiki/Supported-volume-formats#two-level-hierarchy-of-folders).  
An example dataset P30MEMXcreMADMEGFR++L is provided [here](https://doi.org/10.5061/dryad.fj6q57400).  
![The example data structure](/image/dataset.png "Data structure") 

### Stitching
Software: [TeraStitcher](https://abria.github.io/TeraStitcher/)  
 - Input:
     - Reference channel of the dataset as [a two-level hierarchy of folders](https://github.com/abria/TeraStitcher/wiki/Supported-volume-formats#two-level-hierarchy-of-folders). [^1] 
 - Output:  
     - XML descriptor (.xml) containing globally optimal tile positions *(required for **Cell stitching**)*.  
     - Low-resolution (~ 25 µm) stitched volume (.tif) of TIFF 3D format from the *Merge* step *(required for **Registration**)*.  
[^1]: Reference channel - the channel will be used as the reference channel during channel alignment and be used for registration.  

### Channel alignment
Tool: [NM-align](https://github.com/yccc12/COMBINe/tree/main/NuMorph-align)
 - Input:
     - A multi-channel dataset structured as [a two-level hierarchy of folders](https://github.com/abria/TeraStitcher/wiki/Supported-volume-formats#two-level-hierarchy-of-folders).
 - Output:  
     - z displacement table (.mat) if align by translation *(used for **Cell detection**)*.  
     - Translation paramater table (.mat) if align by translation *(used for **Cell detection**)*.  
     - Aligned images (optional for alignment by translation) [^2].  
[^2]: Required if aligned by Elastix.  

### Cell detection
Package: [keras-retinanet](https://github.com/yccc12/COMBINe/tree/main/keras-retinanet)
 - Input:
     - A two-channel dataset structured as [a two-level hierarchy of folders](https://github.com/abria/TeraStitcher/wiki/Supported-volume-formats#two-level-hierarchy-of-folders).  
     - (Optional) z displacement table, translation paramater table, aligned images *(from **Channel alignment**)*. [^3]  
 - Output:  
     - List of predictions (.csv) per tile *(required for **Cell stitching**)*.  
     - Images with predictd boxes (.png).  
[^3]: Channel alignment is optinoal but recommended.  

### Cell stitching
Package: [keras-retinanet](https://github.com/yccc12/COMBINe/tree/main/keras-retinanet)
 - Input:
     - List of predictions (.csv) per tile *(from **Cell detection**)*.  
     - XML descriptor (.xml) containing globally optimal tile positions *(from **Stitching**)*.  
 - Output:  
     - Stitched cell centroids (.csv) per class *(used for **Registration**)*.  

### Registration
Package: [ClearMap2](https://github.com/yccc12/COMBINe/tree/main/ClearMap2)
 - Input:
     - Low-resolution (~ 25 µm) stitched volume (.tif) *(from **Stitching**)*.  
     - Stitched cell centroids (.csv) per class *(from **Cell stitching**)*.   
 - Output:  
     - Annotated cell centroids (.csv) per class *(used for **Analysis**)*.  
     - Registered annotation volume (.mhd) *(used for **Analysis**)*.
     - Voxel-wise cell density counts (.tif) per class *(used for **Analysis**)*.  
     
### Analysis
Script: [Analysis and visualization](https://github.com/yccc12/COMBINe/tree/main/NuMorph-align/analysis)  
 - Input:
     - Annotated cell centroids (.csv) per class *(from **Registration**)*.  
     - Registered annotation volume *(from **Registration**)*.
     - Voxel-wise cell density per class *(from **Registration**)*.  
 - Output:  
     - Table of quantifications per group.  
     - Table of region-wise analysis results.  
     - Voxel-wise analysis results.  

# Contact
ycai23@ncsu.edu; rene.cai@unc.edu

# References
- Bria, Alessandro, and Giulio Iannello. "TeraStitcher-a tool for fast automatic 3D-stitching of teravoxel-sized microscopy images." *BMC bioinformatics* 13.1 (2012): 1-15.  
- Krupa, Oleh, et al. "NuMorph: Tools for cortical cellular phenotyping in tissue-cleared whole-brain images." *Cell reports* 37.2 (2021): 109802.  
- https://github.com/fizyr/keras-retinanet  
- Renier, Nicolas, et al. "Mapping of brain activity by automated volume analysis of immediate early genes." *Cell* 165.7 (2016): 1789-1802.  
