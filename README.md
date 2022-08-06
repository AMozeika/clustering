# clustering
This repository contains the clustering c code,  r code to process results of clustering,  example of data and example of output.  

Here is the instruction (first lines of *.c file) on how to compile and run on multiprocessor machine: 

/*Model based clustering program assuming Gaussian data and using population  dynamics algorithm from the https://arxiv.org/abs/1810.02627 article. This is more efficient paralell version. To transform data uncomment code which calls  the transform_data() function.  To compile on multiprocessor try:  gcc -Wall -fopenmp PopulDynamClustV5v4.c -lm -O3 -o populdynam;  To run on Linux machine use: date; ./populdynam<parameters.in>parameters.out; date; */

The data has to be in the tab separated csv format (and transformed if this is needed). 10dL2c0.csv  is the sample of correlated data (with the same mean vectors and different random covariance matrices) in 10 dimensions and with  2 clusters.   10dL2c0K10.in is the *.in  file for this data. The meaning of the numbers

99191 20000 10 1 10 100 1000 10dL2c0.csv 0 nofile


 in this file is on the first line of the 10dL2c0K10.out file: 

#seed=99191 N=20000   d=10   K1=1   K2=10 rest.=100 t_max=1000   data-file=10dL2c0.csv
#t M1 M2 F
0 9904 10096 24.651864
1 10387 9613 24.154138

So the first number is seed for random number generator, second is the sample size (20000  in this case), third is dimension (d=10), then the range for number of clusters K (K1=1   K2=10 ), number of restarts for the algorithm (here it is rest.=100 but initially use 10), maximum algorithm "time" (here it is t_max=1000 and don't change it for now) and  finally the data-file name (data-file=10dL2c0.csv) . The last two    (0 nofile) are always the same. 

The clustering algorithm produces files which then can be processed  by the r code. For 10dL2c0.csv the r code is 10dL2c0.csv.cluster.statistics.r and 10dL2c0.csv.clustering.statistics.r

The r code takes results of clustering and produces bunch of *.tex and image files for presentations.  Then one need to open 10dL2c0.csv.cluster.statistics.tex and  10dL2c0.csv.clustering.statistics.tex files produced by r code and to compile in tex editor.  This gives pdf presentations  10dL2c0.csv.cluster.statistics.pdf and 10dL2c0.csv.clustering.statistics.pdf

For your data you need to change first lines of r code (see comments in r code). 

Perhaps first you can try to reproduce results of clustering 10dL2c0.csv
