# Efficiently Answering Reliable Shortest Path Queries in Stochastic Road Networks
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

./index [network name] [variance file id] [covariance file id] [query file name]

./ERSP [network name] [variance file id] [covariance file id] [query file name]

./SDRSP [network name] [variance file id] [covariance file id] [query file name]

./SMOGA [network name] [variance file id] [covariance file id] [query file name]

./Maintenance [network name] [variance file id] [covariance file id] [update file name]

./IndIndex [network name] [variance file id] [query file name]

A sample run is ./index NY 3 3, where the program will run all the 10 query files if the last parameter is not specified. The network name could be NY, BAY, or COL. The variance file id can be "1"-"5" as explained in the data description. The covariance id can be "1"-"5" or "K1"-"K4". The query file name can be "dis1"-"dis5" and "alpha1"-"alpha5". The output of the screen will show many auxiliary details, including the index time and memory costs. The results of the query time, the number of hoplinks and path concatenations, will be put in the file with the suffix "Results" in the folder "data/[network name]". The query answers will be in the file named with the suffix "ans".

### Data Description

We stored all the data used in the experiments in the folder data/. The folder contains three sub-folders corresponding to the three networks.

In the folder of each network, the files with suffixes "USA-road-t.[network name].gr", "USA-road-var.[network name].gr", "USA-road-cov.[network name].gr", and "USA-road-d.[network name].gr" corresponds to the mean, variance, covariance, and distance of each edge, respectively. The descriptions can be found in the DIMACS challenge (http://www.diag.uniroma1.it//challenge9/download.shtml). The variance and covariance files with their prefixes "1"-"5" are used in the experiments of varying the CV (as explained in the paper) from 0.1-0.9. The covariance files with their prefixes "K1"-"K4" are used in the experiments of varying K. Note that the variance and covariance files with the prefix "3" use the default setting where K=5 and CV=0.5.

There are also 10 query set files named "dis1"-"dis5" and "alpha1"-"alpha5" corresponding to those stated in the paper. Each query set file contains 1000 lines. Each line contains three values of a CSP query which are s, t, and $\alpha$.

The files with their prefixes "mu" and "sigma" are for the experiments of changing the mean and variance during index maintenance. Each line contains a three values representing the edge id, the multiplier for the mean, and the multiplier for the variance.

The file named "order.txt" is used to build the tree decomposition. One can also generate it by using the function "genorder()" in the index.cpp.