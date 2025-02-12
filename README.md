# scConGraph: a scalable cross-time Context Graph model
scConGraph is a scalable bi-layer graph model that efficiently **integrates cross-time context information**, enabling the comprehensive analysis of **tumor cell dynamic responses** from **paired perturbed or time-seiries single-cell transcriptomics**.

<p align="center">
  <img width="800"  src="Analysis/images/flowchart.jpg">
</p>

## System requirements
scConGraph is currently **available for Linux systems**, as it relies on **LINE**, a `C++`-based embedding method that depends on the GSL package for Linux. If you use the downstream analysis functions in `scConGraph`, which are implemented in Python, there are no system restrictions.

Tested System Configuration:
``` yaml
OS: Linux 3.10.0-1160.el7.x86_64
Python Version: 3.9.19
Processor: x86_64
CPU Cores: 36
Logical CPUs: 36
Total RAM (GB): 251.38
```

## Installation
To ensure a clean and isolated environment for `scConGraph`, we recommend creating a new **Conda environment**:
``` bash
conda create -n scConGraph_env python=3.9 -y
conda activate scConGraph_env
``` 
You can install `scConGraph` directly from PyPI:
``` bash
pip install --index-url https://pypi.org/simple/ scConGraph==1.0.0
```
(Optional) We recommend using `scConGraph` in Jupyter Notebook. To do so, install `ipykernel` and add the environment to Jupyter:
``` bash
pip install ipykernel
python -m ipykernel install --user --name=scConGraph_env --display-name "scConGraph_env"
```
Now, you can select "scConGraph_env" as the kernel when running Jupyter Notebook. 

### Required Dependencies
Importantly, **the LINE toolkit** (LINUX version) must be downloaded and installed from LINE GitHub Repository (https://github.com/tangjianpku/LINE.git) **before using scConGraph**. 



## **Tutorial**
The vignette of `scConGraph` can be found in the project [Wiki](https://github.com/Li-Xinqi/scConGraph/wiki).

## Codes for PDAC Drug Response Analysis
The R scripts used for analyzing drug responses are located in the Analysis/Codes folder. The raw and intermediate data are stored in the Analysis/Data folder. Larger datasets are saved on the [cloud](https://cloud.tsinghua.edu.cn/d/63ff224544874971b0dd/).  
If you have any questions about the codes, please contact the [author](lxq19@mails.tsinghua.edu.cn).
