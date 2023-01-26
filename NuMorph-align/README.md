[NuMorph](https://github.com/ok37/numorph) (MATLAB) is adopted for channel alignment.   

## Installation
This directory contains a light-weight version of NuMorph sufficient for channel alignment. Check [here](https://github.com/ok37/numorph) for the entire toolset.  
Download the codes and run `NM_setup('light')`. If Elastix is needed, run `NM_setup`. 

## Channel alignment

 - Input:
     - A multi-channel dataset structured as [a two-level hierarchy of folders](https://github.com/abria/TeraStitcher/wiki/Supported-volume-formats#two-level-hierarchy-of-folders).
 - Output:  
     - z displacement table (.mat) if align by translation *(used for **Cell detection**)*.  
     - Translation paramater table (.mat) if align by translation *(used for **Cell detection**)*.  
     - Aligned images (optional for alignment by translation) [^1].  
[^1]: Required if aligned by Elastix. 

Add sample information to `./templates/NM_samples.m`.  
Adjust parameters in `./templates/NMp_template`.  
For alignment of a single dataset, follow `NM_run.m`.  
For batch processing, follow `NM_batch.m`.  

## Analysis

 - Input:
     - Annotated cell centroids (.csv) per class *(from **Registration**)*.  
     - Registered annotation volume *(from **Registration**)*.
     - Voxel-wise cell density per class *(from **Registration**)*.  
 - Output:  
     - Table of quantifications per group.  
     - Table of region-wise analysis results.  
     - Voxel-wise analysis results.  
     
Region-wise analysis, follow `./analysis/regionAnaly.m`.  
Voxel-wise analysis, follow `./analysis/voxelAnaly.m`.  
