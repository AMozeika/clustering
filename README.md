# clustering
This repository contains the clustering c code,  r code to process results of clustering,  example of data and example of output.  

The c code implements  paralellised  version of efficient population  dynamics algorithm,  developed for  model based Bayesian  clustering https://arxiv.org/abs/1810.02627,  which  assumes  Gaussian distribution of data.   To transform data uncomment the code which calls  the transform_data() function.  To  compile on the multiprocessor (Linux) machine the command  "gcc -Wall -fopenmp PopulDynamClustV5v4.c -lm -O3 -o populdynam;" was used in the terminal.  To run on the multiprocessor (Linux) machine the command  "date; ./populdynam<parameters.in>parameters.out; date;" was used in the terminal.  

The data has to be in the tab separated csv format (and transformed if this is needed). 10dL2c0.csv  is the sample of correlated data (with the same mean vectors and different random covariance matrices) in 10 dimensions and with  2 clusters.   10dL2c0K10.in is the *.in  file for this data. The meaning of the numbers "99191 20000 10 1 10 100 1000 10dL2c0.csv 0 nofile" in this file is given  in  the first line of the 10dL2c0K10.out file: 

#seed=99191 N=20000   d=10   K1=1   K2=10 rest.=100 t_max=1000   data-file=10dL2c0.csv
#t M1 M2 F
0 9904 10096 24.651864
1 10387 9613 24.154138

So first number is the seed for random number generator, second is the sample size "N=20000", third is dimension (d=10), then the range for number of clusters to consider "K1=1   K2=10", number of restarts for the algorithm (here it is rest.=100), maximum algorithm "time" (here it is t_max=1000) and  finally the data-file name (data-file=10dL2c0.csv) . The last two    "0 nofile" are always the same. 

The clustering algorithm produces files which then can be processed  by the r code. For 10dL2c0.csv the r code is 10dL2c0.csv.cluster.statistics.r and 10dL2c0.csv.clustering.statistics.r

The r code takes results of clustering and produces bunch of *.tex and image files for presentations.  Then one need to open 10dL2c0.csv.cluster.statistics.tex and  10dL2c0.csv.clustering.statistics.tex files produced by r code and to compile in tex editor.  This gives pdf presentations  10dL2c0.csv.cluster.statistics.pdf and 10dL2c0.csv.clustering.statistics.pdf

For your data you need to change first lines of r code (see comments in r code). 
