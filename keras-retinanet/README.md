[fizyr/keras-retinanet](https://github.com/fizyr/keras-retinanet) (Python) is adopted for cell detection.  

## Installation with Anaconda (Windows)

Download the codes from [fizyr/keras-retinanet](https://github.com/fizyr/keras-retinanet).   
Replace `./keras_retinanet/utils/image.py`, `./keras_retinanet/utils/colors.py`, and `./keras_retinanet/utils/gpu.py`.  

Check [here](https://docs.anaconda.com/anaconda/install/index.html) for Anaconda documentation.  
Open Anaconda Prompt and create an virtual environment - the package versions below are not mandontary but should work.  
`conda create -n keras-retinanet python=3.7 tensorflow-gpu=2.3 tensorflow=2.3=mkl_py37h936c3e2_0 keras=2.4 numpy=1.19.2`      
Activate the environment: `conda activate keras-retinanet`.  

How I had to do it:
```bash
conda create -n keras-retinanet python=3.7
conda activate keras-retinanet
pip install tensorflow==2.3
```

Go to the code directory, e.g. `.../keras-retinanet-main/`, in Anaconda Prompt.  
Install keras-retinanet `pip install . --user`.  
To run the code directly from the directory, run `python setup.py build_ext --inplace` to compile Cython code.

## Training and testing
Please refer to [fizyr/keras-retinanet](https://github.com/fizyr/keras-retinanet) for more details.  

## Cell detection

 - Input:
     - A two-channel dataset structured as [a two-level hierarchy of folders](https://github.com/abria/TeraStitcher/wiki/Supported-volume-formats#two-level-hierarchy-of-folders).  
     - (Optional) z displacement table, translation paramater table, aligned images *(from **Channel alignment**)*. [^1]  
 - Output:  
     - List of predictions (.csv) per tile *(required for **Cell stitching**)*.  
     - Images with predictd boxes (.png).  
[^1]: Channel alignment is optinoal but recommended. 

The model is trained to detect six classes of cells in MADM mouse brains, red, green, yellow: neuron and astrocyte.  

## Cell stitching

 - Input:
     - List of predictions (.csv) per tile *(from **Cell detection**)*.  
     - XML descriptor (.xml) containing globally optimal tile positions *(from **Stitching**)*.  
 - Output:  
     - Stitched cell centroids (.csv) per class *(used for **Registration**)*.  

Follow `inferenceStitch.py`, an automated cell detection and stitching pipeline designed for a two-channel [TeraStitcher compatible](https://github.com/abria/TeraStitcher/wiki/Supported-volume-formats#two-level-hierarchy-of-folders) dataset.
