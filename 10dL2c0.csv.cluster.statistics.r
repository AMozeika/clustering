data_file<-'10dL2c0.csv';  #insert name of data file here
#
label<-c('x1',	'x2', 	'x3', 	'x4', 	'x5', 	'x6',  	'x7', 	'x8', 'x9', 	'x10'); #labels for columns in data
#
d=10; #dimension of data
#
K1=1;  K2=10;#range for the number of clusters
#
COLOR=1; #1 for color histograms and 0 for black/grey histograms
#
Color<-rainbow(K2); #array of colors
#
library(dendextend);
#
file<-paste(data_file,'.K',K1,'.clusters.csv', sep='');
#
data<-read.table(file,  header = FALSE, sep = "\t"); head(data); #import clustering result
#
mat <- matrix(, nrow = 2, ncol = d); #define matrix to store max and min of columns
# 
for (j in 1:d){          mat[1,j]=max(data[,j]); mat[2,j]=min(data[,j])     };   #compute max and min of columns
#
slides_file<-paste(data_file, '.cluster.statistics.tex',sep='');
cat('\\documentclass{beamer}\n ',file=slides_file, sep=' ');
cat('\\usepackage{graphicx}\n ',file=slides_file, sep=' ', append=TRUE);
cat('\\mode<presentation>\n',file=slides_file, sep=' ', append=TRUE);
cat('\\usetheme{Warsaw}\n',file=slides_file, sep=' ', append=TRUE);
cat('\\title[', data_file, ' ]{Cluster  Statistics} \n ',file=slides_file, sep=' ', append=TRUE);
cat('\\author{Population dynamics V 5.4} \n ',file=slides_file, sep=' ', append=TRUE);
cat('\\date{\\today}  \n ',file=slides_file, sep=' ', append=TRUE);
cat('\\setbeamertemplate{page number in head/foot}[totalframenumber]\n\\begin{document} \n ',file=slides_file, sep=' ', append=TRUE);
cat('\\begin{frame}\\titlepage\\end{frame}\n ',file=slides_file, sep=' ', append=TRUE);
cat('%list of slides\n',file=slides_file, sep=' ', append=TRUE);
#
for (K in K1:K2)
{
	file<-paste(data_file,'.K',K,'.clusters.csv', sep='');
	clusters<-read.table(file,  header = FALSE, sep = "\t"); head(clusters); #import clustering result
	X <- split(clusters, clusters[d+1]);  #use values in column d+1 to split the data into K data-sets
	Y <- lapply(seq_along(X), function(x) as.data.frame(X[[x]])[, 1:d]);  #remove d+1 column
	cluster <- list(); #create list to store clusters
	for (i in 1:K){ cluster[[i]] <- Y[[i]];} #populate list
	#
	#
	####################begin create slide with 2d projection of K clusters#################
	x=1; y=2; #choose two columns from data
	slide<-paste(data_file,'.',label[x],'.',label[y],'.K', K,'.tex', sep='');
	sink(slide);
	figure<-paste(label[y],'vs',label[x],'K',  K, '.png', sep='');
	png(figure);
	plot(cluster[[1]][,x],cluster[[1]][,y],pch=1,col=1, xlab=label[x], ylab=label[y], xlim=c(  min(   data[x]  ),  max(  data[x]   )  ),  ylim=c(  min(   data[y]   ),  max(   data[y]    )    )     );
	if(K>1)
	{
		for (i in 2:K)
		{ 
			par(new=TRUE);
			if(i<8)
			{
				plot(cluster[[i]][,x],cluster[[i]][,y],pch=i,col=i, xlab='', ylab='',xlim=c(  min(   data[x]      ),  max(   data[x]     )  ),  ylim=c(  min(    data[y]      ),  max(    data[y]      )  ) ,axes=F)	
			}
			else
			{
				plot(cluster[[i]][,x],cluster[[i]][,y],pch=i,col=Color[i], xlab='', ylab='',xlim=c(  min(   data[x]      ),  max(   data[x]     )  ),  ylim=c(  min(    data[y]      ),  max(    data[y]      )  ) ,axes=F)	
			}		
		}
	}
	dev.off();
	#
	width=110; height=60; #figure dimensions
	#
	cat('\\begin{frame}\n','\\frametitle{Results for $K=', K, '$.}\n\\begin{figure}[h]\\setlength{\\unitlength}{1mm}\n\\begin{picture}(',width,',', height,')\n',file=slide, sep=' ',append=TRUE);
	#
	width=height;  height=height; #picture dimensions
	#
	cat('\\put(0,0){\\includegraphics[height=',height,'\\unitlength, width=',width,'\\unitlength]{',figure, '}}\n',file=slide, sep='',append=TRUE);
	cat('\\end{picture}\n \\end{figure}\n\\end{frame}\n', file=slide, sep=' ',append=TRUE);
	sink();
	cat('\\input{',slide, '}\n', file=slides_file, sep='', append=TRUE);
	#####################end create slide with 2d projection of K clusters###################
	#
	#########################begin create dendrogram#########################
	#
	if(K>2)
	{
		#compute matrix of K-L distances for clusters K
		D <- matrix(0, nrow = K, ncol = K);
		for (r in 1:K)
		{
			for (c in 1:K) 
			{ 	
				if(r!=c)
				{	
					#compute covariances	of clusters 1 and 2 
					Cov1=var(cluster[[r]]);
					Cov2=var(cluster[[c]]);	
					#compute means of clusters 1 and 2
					m1<-c(); 
					m2<-c(); 
					for (i in 1:d)
					{
						m1[i]=mean(cluster[[r]][,i])
						m2[i]=mean(cluster[[c]][,i])
					};
					#comp. inverse of Cov2 
					Pr2<-solve(Cov2); 
					#compute quadr. form
					QF=( t(m2-m1)%*%Pr2%*%(m2-m1));
					#use above to compute
					D[r,c]=0.5*(sum(diag(Pr2%*%Cov1)) + QF[1,1] + log(det(Cov2)/det(Cov1))-d);
				}
			}  
		};
	#
	#use K-L distance to compute J-S distance
	JD <- matrix(0, nrow = K, ncol = K);
	for (r in 1:(K-1))
	{
		for (c in (r+1):K)
		{ 	
			JD[r,c]=0.5*(D[r,c]+D[c,r]);
			JD[c,r]=JD[r,c];
			    
		}  
	};
	#
	Dist <- as.dist(JD, diag = TRUE);
	hclusters <- hclust(Dist);
	dend <- as.dendrogram(hclusters);
	if(K<8)
	{
		colors_to_use <-1:K;
	}
	else
	{
		colors_to_use <- 1:K;
		for (i in 8:K)
		{
			colors_to_use[i] <- Color[i];
		}
	}
	colors_to_use <- colors_to_use[order.dendrogram(dend)];
	labels_colors(dend) <-colors_to_use;
	figure<-paste('dendr','K',  K, '.png', sep='');
	png(figure);
	#plot(hclusters);	
	plot(dend);
	dev.off();	
	#create slide 
	slide<-paste(data_file, '.K', K,'.dendr.tex', sep='');
	sink(slide);
	width=110; height=70; #figure dimensions
	#
	cat('\\begin{frame}\n','\\frametitle{Results for $K=', K, '$. Kullback--Leibler divergence.}\n\\begin{figure}[h]\\setlength{\\unitlength}{1mm}\n\\begin{picture}(',width,',', height,')\n',file=slide, sep=' ',append=TRUE);
	#
	width=height;  height=height; #picture dimensions
	#
	cat('\\put(0,0){\\includegraphics[height=',height,'\\unitlength, width=',width,'\\unitlength]{',figure, '}}\n',file=slide, sep='',append=TRUE);
	cat('\\end{picture}\n \\end{figure}\n\\end{frame}\n', file=slide, sep=' ',append=TRUE);
	sink();
	cat('\\input{',slide, '}\n', file=slides_file, sep='', append=TRUE);
}
#
if(K>2)
{
	Clusters<-hclusters$order;	
}	
else
{
	Clusters<-1:K;		
}
#####################################end create dendrogram#####################################################################
#
for (i in 1:K)
{
	#generate figures of histograms of marginals
	slide<-paste(data_file, '.K', K,'.cluster.', Clusters[i], '.tex', sep='');
	sink(slide);
	width=110; height=60; #figure dimensions in mm
	cat('\\begin{frame}\n','\\frametitle{Results for $K=', K, '$. Cluster $',Clusters[i],'$.}\n\\begin{figure}[h]\\setlength{\\unitlength}{1mm}\n\\begin{picture}(',width,',', height,')\n',file=slide, sep=' ',append=TRUE);
	#		
	nr=2; nc=5;		 #nr*nc=d
	#histogram dimensions
	hy=(height/nr)*0.5;
	hx=hy; 
	#
	j=1;
	for (r in (nr-1):0)
	{
		for (c in 0:(nc-1))
		{
			figure<-paste('K',  K,    'clust', Clusters[i], label[j],'.png', sep='');
			png(figure);
			x<-cluster[[Clusters[i]]][,j];
			if(Clusters[i]<8)
			{
				if (COLOR>0)
				{
					hist(data[,j], prob=TRUE, xlab=label[j], ylab='', main='',col=8, xlim=c( mat[2,j]  ,  mat[1,j] ) ,cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3); 
					hist(x, prob=TRUE, col=Clusters[i],  add=T);	
				}
				else
				{
					hist(data[,j], prob=TRUE, xlab=label[j], ylab='', main='',col=8, xlim=c( mat[2,j]  ,  mat[1,j] ) ,cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3); 
					hist(x, prob=TRUE, col=1,  add=T);							
				}
			}
			else
			{
						if (COLOR>0)
						{
							hist(data[,j], prob=TRUE, xlab=label[j], ylab='', main='',col=8, xlim=c( mat[2,j]  ,  mat[1,j] ) ,cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3); 
							hist(x, prob=TRUE, col=Color[Clusters[i]],  add=T);	
						}
						else
						{
							hist(data[,j], prob=TRUE, xlab=label[j], ylab='', main='',col=8, xlim=c( mat[2,j]  ,  mat[1,j] ) ,cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3); 
							hist(x, prob=TRUE, col=1,  add=T);							
						}
			}
			dev.off();
			space=4; #blank space between figures
			cat('\\put(',c*(hx+space), ',' , r*(hy+space), '){\\includegraphics[height=',hy,'\\unitlength, width=',hx,'\\unitlength]{',figure, '}}\n',file=slide, sep='',append=TRUE);
			j=j+1;
		}
	}			
			#
			cluster_size<-paste('cluster size=',nrow(cluster[[Clusters[i]]]),';', sep='');
			text<-paste(cluster_size, sep=' ');
			cat('\\end{picture}\n \\begin{center}\n', '\\small{',text, '}','\\end{center}\n\\end{figure}\n\\end{frame}\n', file=slide, sep=' ',append=TRUE);
			sink();
			#
			cat('\\input{',slide, '}\n', file=slides_file, sep='', append=TRUE);
			#
#			
	}
}
cat('\\end{document} \n ',file=slides_file, sep=' ', append=TRUE);
#
cat('Open', slides_file, 'file\n', sep=' ');

