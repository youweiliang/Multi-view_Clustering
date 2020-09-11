## Multi-view Spectral Clustering Algorithms
This repository contains 7 multi-view spectral clustering algorithms (and a single-view spectral clustering algorithm) used for comparison in our paper  "[Multi-view Graph Learning by Joint Modeling of Consistency and Inconsistency](https://arxiv.org/abs/2008.10208)", which is the follow-up work of our ICDM paper "Consistency Meets Inconsistency: A Unified Graph Learning Framework for Multi-view Clustering". 
The code of some algorithms was gathered from the websites of the authors of the original papers and was later fixed and optimized by us. 
Please see our paper for the details of these algorithms (the folder names correspond to the abbreviations for the algorithms in our paper).


### Preparation
* **Windows 64bit**: 
Add some helper files to MATLAB path by `addpath('MinMaxSelection')` command in MATLAB command window.
* **Linux, Windows 32bit and Mac OS**: 
Add some helper files to MATLAB path by `addpath('MinMaxSelection')` command in MATLAB command window. Then recompile the helper functions by running `minmax_install`.


### Example usage
The file `test.m` contains examples to use all the algorithms. 
