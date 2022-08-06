data_file<-'10dL2c0.csv';  #name of data file
K1=1;  K2=10;#number of clusters
slides_file<-paste(data_file, '.clustering.statistics.tex',sep='');
cat('\\documentclass{beamer}\n ',file=slides_file, sep=' ');
cat('\\mode<presentation>\n',file=slides_file, sep=' ', append=TRUE);
cat('\\usetheme{Warsaw}\n',file=slides_file, sep=' ', append=TRUE);
cat('\\title[', data_file, ' ]{Clustering Statistics} \n ',file=slides_file, sep=' ', append=TRUE);
cat('\\author{Population dynamics V 5.4} \n ',file=slides_file, sep=' ', append=TRUE);
cat('\\date{\\today}  \n ',file=slides_file, sep=' ', append=TRUE);
cat('\\setbeamertemplate{page number in head/foot}[totalframenumber]\n\\begin{document} \n ',file=slides_file, sep=' ', append=TRUE);
cat('\\begin{document} \n ',file=slides_file, sep=' ', append=TRUE);
cat('\\begin{frame}\\titlepage\\end{frame}\n ',file=slides_file, sep=' ', append=TRUE);
#
F_min <- matrix(0, nrow = K2, ncol = 2); #define matrix to store F min  and K
#
for (K in K1:K2)
{
	#import data from file	
	file<-paste(data_file ,'.K',K,'.F.csv', sep='');
	data<-read.table(file,  header = T, sep = "\t"); head(data); #import data 
	F_min[K,1]=K;
	F_min[K,2]=min(data$F);
}
#
#####################begin generate F vs K figure and slide######################################
figure<-paste('FvsK.png', sep='');
png(figure);
plot(F_min[,1], F_min[,2], xlab='K', ylab='min F',pch=3, panel.first = grid(lwd=2));
dev.off();
#
#generate slide 
file<-paste('FvsK.tex', sep='');
cat('\\input{',file, '}\n', file=slides_file, sep='', append=TRUE);
#
sink(file);
cat('\\begin{frame}\\frametitle{ $\\mathrm{F}$ as a function of $K$}\\begin{center}\\includegraphics[height=7cm]{',figure,'} \\end{center}\\end{frame}',sep='');
sink();
#####################end generate F vs K figure and slide#############################################################################
#
####################begin generate F vs ln(K) figure and slide###########################################################################
figure<-paste('FvslnK.png', sep='');
png(figure);
plot(F_min[,1], F_min[,2]+log(F_min[,1]), xlab='K', ylab='min F+ln(K)',pch=3, panel.first = grid(lwd=2));
dev.off();
#
#generate slide 
file<-paste('FvslnK.tex', sep='');
cat('\\input{',file, '}\n', file=slides_file, sep='', append=TRUE);
#
sink(file);
cat('\\begin{frame}\\frametitle{ $\\mathrm{F}+\\ln (K)$ as a function of $K$}\\begin{center}\\includegraphics[height=7cm]{',figure,'} \\end{center}\\end{frame}',sep='');
sink();
####################end generate F vs ln(K) figure and slide###########################################################################
#
####################begin generate dF vs K figure and slide###########################################################################
dF <- matrix(0, nrow = (K2-1), ncol = 2); #define matrix to store dF min  and K
for (K in 1:(K2-1)){dF[K,2]=-1*(F_min[K+1,2]+log(K+1)-F_min[K,2]-log(K));      dF[K,1]=K+1;  };
figure<-paste('dFvsK.png', sep='');
png(figure);
plot(dF[,1], dF[,2], xlim=c(1, K2), xlab='K', ylab='dF',pch=3, panel.first = grid(lwd=2));
dev.off();
#
#generate slide 
file<-paste('dFvsK.tex', sep='');
cat('\\input{',file, '}\n', file=slides_file, sep='', append=TRUE);
#
sink(file);
cat('\\begin{frame}\\frametitle{ $\\mathrm{d}\\mathrm{F}$ as a function of $K$}\\begin{center}\\includegraphics[height=7cm]{',figure,'} \\end{center}\\end{frame}',sep='');
sink();
####################end generate dF vs K figure and slide###########################################################################
#
for (K in K1:K2)
{
	if(K>1)	
	{
		#import data from file	
		file<-paste(data_file ,'.K',K,'.F.csv', sep='');
		data<-read.table(file,  header = T, sep = "\t"); head(data); #import data 
		#
		#generate histogram
		figure<-paste('K',  K,    'F', '.png', sep='');
		png(figure);
		hist(data$F, prob=FALSE, xlab='F', main='');
		dev.off();
		#
		#generate slide 
		file<-paste('K', K,'F', '.tex', sep='');
		cat('\\input{',file, '}\n', file=slides_file, sep='', append=TRUE);
		#
		minF<-paste("$\\small{\\min{ \\mathrm{F}}}=", min(data$F), '$',sep='');
		#
		sink(file);
		cat('\\begin{frame}\\frametitle{Histogram of F for K=',K ,'}\\begin{center}\\includegraphics[height=7cm]{',figure,'} \\end{center}',minF,'\\end{frame}',sep='');
		sink();
	}
}
cat('\\end{document} \n ',file=slides_file, sep=' ', append=TRUE);
#
cat('Open', slides_file, 'file\n', sep=' ');

