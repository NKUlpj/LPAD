# LPAD
Detection of Topological Associated Domains from Hi-C Data using Network Construction and Label propagation

## 1. Hi-C Data 

In our experiment, we used 3 three normalized Hi-C matrix pre-processed by Bing Ren's Lab. Data is available at : http://chromosome.sdsc.edu/mouse/hi-c/download.html.
Besides,we also used a simulated data extract from chr6 on mouse cortex cell line, normalized Matrix is in the folder  `data`.

## 2.Usage:

run the PLA_2.py

Parameters are as follows:

 * filepath: the path of a `intra-chromosomal` Hi-C matrix seperated by Tab with N by N shaped
 * w: window-size,the default is 5.
 * c: the restart probability, the default is 0.9
 * outpath: the storage path and filename of result, the default is ./out
 * p (<b>optional</b>): the Penalty coefficient, the default is 1.

### Requirements
 * numpy
 * matplotlib
 * seaborn
 * colormap
 * pyrwr: we use this package download from https://github.com/jinhongjung/pyrwr and modify the source code. Therefore, we uploaded the code in this repository.
 * scipy
 * tdqm
 * pandas

 ## 3.Output

LPAD produces a files. Note that,The output file will be named by the user in which list the TADs extracted from the input data. The first column is the start bin, and the second column is the stop bin.
