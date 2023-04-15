# Efficiently Answering Reliable Shortest Path Queries in Stochastic Road Networks
========================================================================
This repository contains the source codes and data for this paper. 

Usage
---------------

### Environment

g++ version: 8.4.0 

OS: Ubuntu/CentOS (Unix)

### Compilation

cd code

make

### Execution

In the folder code/,

./index [network name] [variance path] [covariance path]

### Data Description

The data used in the experiments are stored in the folder data/. The folder contains three sub-folders corresponding to the three networks.

In the folder of each network, it contains two files named "USA-road-t.[network name].gr" and "USA-road-d.[network name].gr" corresponding to the travel time and the distance of each edge. The descriptions can be found in the DIMACS challenge (http://www.diag.uniroma1.it//challenge9/download.shtml).

It also contains 10 query set files named dis1-dis5 and alpha1-alpha5 corresponding to those stated in the paper. Each query set file contains 1000 lines. Each line contains three values of a CSP query which are s, t, and C.

The file named "order.txt" is used to build the tree decomposition. One can also generate it by using the function "genorder()" in the index.cpp.
