# scConGraph: a scalable cross-time Context Graph model
scConGraph is a scalable bi-layer graph model that efficiently **integrates cross-time context information**, enabling the comprehensive analysis of **tumor cell dynamic responses** from **paired perturbed or time-seiries single-cell transcriptomics**.

<p align="center">
  <img width="800"  src="https://github.com/Li-Xinqi/scConGraph/assets/53567070/bf948041-ed83-4df8-b487-ebe81c6e9a43">
</p>

## System requirements
scConGraph is now available for Linux, with a Windows version planned for recent release. scConGraph requires only a standard computer with enough RAM to perform in-memory computations.

## Installation
CINEMA-OT requires `python` version 3.7+.  Install directly from pip with:

    pip install scConGraph
    
The LINE toolkit achieved in LINUX version should be downloaded and installed from https://github.com/tangjianpku/LINE.git before using scConGraph. 

## Tutorial
The vignette of `scConGraph` can be found in the project [wiki](https://github.com/Li-Xinqi/scConGraph/wiki).

## Codes for Drug Response Analysis
The R scripts used for analyzing drug responses are located in the Analysis/Codes folder. The raw and intermediate data are stored in the Analysis/Data folder. Larger datasets are saved on the [cloud](https://cloud.tsinghua.edu.cn/d/63ff224544874971b0dd/).  
If you have any questions about the codes, please contact the author at lxq19@mails.tsinghua.edu.cn.
