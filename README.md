## Multi-view Spectral Clustering Algorithms
**This repository contains MATLAB code for 7 multi-view spectral clustering algorithms** (and a single-view spectral clustering algorithm) used for comparison in our ICDM paper "[Consistency Meets Inconsistency: A Unified Graph Learning Framework for Multi-view Clustering](https://youweiliang.github.io/publications/icdm19)". 
The code of some algorithms was gathered from the websites of the authors of the original papers and was later fixed and optimized by us. 
Please see our paper for the details of these algorithms (the folder names correspond to the abbreviations for the algorithms in our paper, namely, **AASC, AWP, CoReg, MCGC, MVGL, RMSC,** and **WMSC**). In each of these folders, there is a main file `xxx_main.m` for the algorithm where `xxx` is the algorithm name.  

The original papers for the 7 multi-view spectral clustering algorithms and the single-view spectral clustering (SC) algorithm are:
* Huang et al., 2012. Affinity Aggregation for Spectral Clustering
* Nie et al., 2018. Multiview Clustering via Adaptively Weighted Procrustes
* Kumar et al., 2011. Co-regularized Multi-view Spectral Clustering
* Zhan et al., 2018. Multiview Consensus Graph Clustering
* Zhan et al., 2017. Graph Learning for Multiview Clustering
* Xia et al., 2014. Robust Multi-view Spectral Clustering via Low-rank and Sparse Decomposition
* Zong et al., 2018. Weighted Multi-view Spectral Clustering Based on Spectral Perturbation
* Ng et al., 2002. On Spectral Clustering: Analysis and an Algorithm


<!-- **For the code of our Multi-view Graph Learning algorithm, please see [this repository](https://github.com/youweiliang/Multi-view_Graph_Learning).** -->

### Dataset
All datasets used in our paper are available at [Baidu Cloud](https://pan.baidu.com/s/1bAfDcgH3NguqWM6saDTv1g) with code `pqti` and [Google Drive](https://drive.google.com/drive/folders/1UtjL0Og7ALs9AJq9XnkdrYUmr5rudCyk?usp=sharing). Each dataset is a mat file containing 2 variables `fea` (i.e., a MATLAB cell of features) and `gt` (i.e., ground truth label), except the file `flower17.mat` which contains a cell of distance matrices and ground truth since features are unavailable. 
* The distance matrices in `flower17.mat` should be squared before passing them into the SGF and DGF functions, and the string `original` should be passed into the functions as the metric parameter. 
* The datasets `Reuters`, `Reuters-21578`, `BBCSport`, and `CiteSeer` are text datasets with word frequence as features and thus should be used with the `cosine` metric for computing distance matrices. 

### Preparation
* **Windows 64bit**: 
Add some helper files to MATLAB path by `addpath('MinMaxSelection'); addpath('utils')` command in MATLAB command window.
* **Linux, Windows 32bit and Mac OS**: 
Add some helper files to MATLAB path by `addpath('MinMaxSelection'); addpath('utils')` command in MATLAB command window. Then recompile the helper functions by running `minmax_install`.


### Example usage
The file `test.m` contains examples to use all the algorithms. 
