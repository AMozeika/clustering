/*Model based clustering program assuming Gaussian data and using population  dynamics algorithm from the https://arxiv.org/abs/1810.02627 article. This is more efficient paralell version. To transform data uncomment code which calls  the transform_data() function.  To compile on multiprocessor try:  gcc -Wall -fopenmp PopulDynamClustV5v4.c -lm -O3 -o populdynam;  To run on Linux machine use: date; ./populdynam<parameters.in>parameters.out; date; */
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <omp.h>

void nrerror(error_text)
char error_text[];
{
void exit();
     fprintf(stderr,"Numerical Recipes run-time error...\n");
     fprintf(stderr,"%s\n",error_text);
     fprintf(stderr,"...now exiting to system...\n");
     getchar();
     exit(1);
}

int *ivector(nl,nh)
int nl,nh;
{
int *v;
     v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
     if (!v) nrerror("allocation failure in ivector()");
     return v-nl;
}

short *svector(nl,nh)
int nl,nh;
{
short *v;
     v=(short *)malloc((unsigned) (nh-nl+1)*sizeof(short));
     if (!v) nrerror("allocation failure in svector()");
     return v-nl;
}

double *dvector(nl,nh)
int nl,nh;
{
double *v;
     v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
     if (!v) nrerror("allocation failure in dvector()");
     return v-nl;
}

double **dmatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
int i;
double **m;
     m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
     if (!m) nrerror("allocation failure 1 in dmatrix()");
     m -= nrl;
     for(i=nrl;i<=nrh;i++) {
          m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
          if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
          m[i] -= ncl;
     }
     return m;
}


int **imatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
int i,**m;
     m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
     if (!m) nrerror("allocation failure 1 in imatrix()");
     m -= nrl;
     for(i=nrl;i<=nrh;i++) {
          m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
          if (!m[i]) nrerror("allocation failure 2 in imatrix()");
          m[i] -= ncl;
     }
     return m;
}

short **smatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
int i;
short **m;
     m=(short **)malloc((unsigned) (nrh-nrl+1)*sizeof(short*));
     if (!m) nrerror("allocation failure 1 in imatrix()");
     m -= nrl;
     for(i=nrl;i<=nrh;i++) {
          m[i]=(short *)malloc((unsigned) (nch-ncl+1)*sizeof(short));
          if (!m[i]) nrerror("allocation failure 2 in smatrix()");
          m[i] -= ncl;
     }
     return m;
} 


void free_ivector(v,nl,nh)
int *v,nl,nh;
{
     free((char*) (v+nl));
}

void free_svector(v,nl,nh)
short *v;
int nl,nh;
{
     free((char*) (v+nl));
}

void free_imatrix(m,nrl,nrh,ncl,nch)
int **m;
int nrl,nrh,ncl,nch;
{
int i;
     for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
     free((char*) (m+nrl));
}

void free_dvector(v,nl,nh)
double *v;
int nl,nh;
{
     free((char*) (v+nl));
}

void free_dmatrix(m,nrl,nrh,ncl,nch)
double **m;
int nrl,nrh,ncl,nch;
{
int i;
     for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
     free((char*) (m+nrl));
}

void free_smatrix(m,nrl,nrh,ncl,nch)
short **m;
int nrl,nrh,ncl,nch;
{
int i;
     for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
     free((char*) (m+nrl));
}


//global variables
int K, d, N; 



int power(int base, unsigned int exp){

    if (exp == 0)
        return 1;
    int temp = power(base, exp/2);
    if (exp%2 == 0)
        return temp*temp;
    else
        return base*temp*temp;

}



void load_data(char *name, double  **x)
{
	FILE *fopen(),*fp;
	int i, j;
	fp=fopen(name,"r");
	
   
   for(i=1; i<=N; i++)
   {
   
  		 for(j=1; j<=d; j++)
   				fscanf(fp,"%lf", &x[i][j]); 
   
   
   }
   
  
    
   fclose(fp);
}


#define TINY 1.0e-20;

double dabs(x) double x; { if(x<0.0) return(-1.0*x); return(x); }

int dludcmp(double **a,int n,int *indx,double *d)
{
int i,imax=0,j,k;
double big,dum,sum,temp,dabs();
double *vv,*dvector();
void nrerror();
void free_dvector();
     vv=dvector(1,n);
     *d=1.0;
     for (i=1;i<=n;i++) {
          big=0.0;
          for (j=1;j<=n;j++)
               if ((temp=dabs(a[i][j])) > big) big=temp;
          if (big == 0.0) {
              puts("ERROR: Singular matrix in routine dLUDCMP");
              /* nrerror("Singular matrix in routine dLUDCMP"); */
              
              //printf("row i=%d\n",i);
              //for(j=1;j<=n;j++) printf("a[%d,%d]=%lf ",i,j,a[i][j]);
              //puts(""); getchar();
              return (-1);
              }
              
          vv[i]=1.0/big;
          }
     for (j=1;j<=n;j++) {
          for (i=1;i<j;i++) {
               sum=a[i][j];
               for (k=1;k<i;k++) sum -= (a[i][k]*a[k][j]);
               a[i][j]=sum;
               }
          big=0.0;
          for (i=j;i<=n;i++) {
               sum=a[i][j];
               for (k=1;k<j;k++)
                    sum -= (a[i][k]*a[k][j]);
               a[i][j]=sum;
               if ( (dum=vv[i]*dabs(sum)) >= big) {
                    big=dum;
                    imax=i;
                  }
               }
          if (j != imax) {
               for (k=1;k<=n;k++) {
                    dum=a[imax][k];
                    a[imax][k]=a[j][k];
                    a[j][k]=dum;
                    }
               *d = -1.0*(*d);
               vv[imax]=vv[j];
               }
          indx[j]=imax;
          if (a[j][j] == 0.0) a[j][j]=TINY;
          if (j != n) {
               dum=1.0/(a[j][j]);
               for (i=j+1;i<=n;i++) a[i][j] *= dum;
               }
          }
     free_dvector(vv,1,n);
	 return(0);
}

void cpy_dmat(double **dest, double **srce, int n)
{
int r,c;
   for(r=1;r<=n;r++) for(c=1;c<=n;c++) dest[r][c]=srce[r][c];
}



//return an error for singular matrix
double det_dmat(double **mat,int n,int *error)
{
double d, **dum, **dmatrix();
int j, *indx, *ivector();
void free_dmatrix(), free_ivector(); int dludcmp(); void cpy_dmat();

     dum = dmatrix(1,n,1,n);
     indx = ivector(1,n);
     cpy_dmat(dum,mat,n);
#if 0
     dludcmp(dum,n,indx,&d);
               /* printf("      d=%lf\n",d); fflush(stdout); */
     for(j = 1; j <= n; j++) {
        d *= (dum[j][j]);
        /* printf("      d=%lf\n",d); fflush(stdout); */
        }
#else
     if (dludcmp(dum,n,indx,&d)==0) {
                   /* printf("      d=%lf\n",d); fflush(stdout); */
         for(j = 1; j <= n; j++) {
            d *= (dum[j][j]);
            /* printf("      d=%lf\n",d); fflush(stdout); */
            }

         *error=0;
     }
     else {
         *error=1;
     }
#endif
     
     free_ivector(indx,1,n);
     free_dmatrix(dum,1,n,1,n);
     return(d);
}

#if 1
//return an error for singular matrix
int inv_dmat(double **inv,double **mat,int n)
{
    double d, **dum, **dmatrix(), *col, *dvector();
    int i, j, error=0, *indx, *ivector();
    void free_dmatrix(), free_dvector(), free_ivector();
    int dludcmp(), dlubksb(); void cpy_dmat();

    dum = dmatrix(1,n,1,n);
    indx = ivector(1,n);
    col = dvector(1,n);
    cpy_dmat(dum,mat,n);

    if (dludcmp(dum,n,indx,&d)==0) {
        for(j = 1; j <= n; j++) {
            for(i = 1; i <= n; i++) col[i] = 0.0;
            col[j] = 1.0;
            dlubksb(dum,n,indx,col);
            for(i = 1; i <= n; i++) inv[i][j] = col[i];
        }
    }
    else {
        error=1;      
    }

    free_dvector(col,1,n);
    free_ivector(indx,1,n);
    free_dmatrix(dum,1,n,1,n);

    return (error);
}
#else
int inv_dmat(double **inv,double **mat,int n)
{
double d, **dum, **dmatrix(), *col, *dvector();
int i, j, *indx, *ivector();
void free_dmatrix(), free_dvector(), free_ivector();
int dludcmp(), dlubksb(), cpy_dmat();
     dum = dmatrix(1,n,1,n);
     indx = ivector(1,n);
     col = dvector(1,n);
     cpy_dmat(dum,mat,n);
     dludcmp(dum,n,indx,&d);
     for(j = 1; j <= n; j++) 
        { for(i = 1; i <= n; i++) col[i] = 0.0;
          col[j] = 1.0;
          dlubksb(dum,n,indx,col);
          for(i = 1; i <= n; i++) inv[i][j] = col[i];
        }
     free_dvector(col,1,n);
     free_ivector(indx,1,n);
     free_dmatrix(dum,1,n,1,n);
   return(0);
}
#endif

int dlubksb(double **a,int n,int *indx,double *b)
{
int i,ii=0,ip,j;
double sum;
     for (i=1;i<=n;i++) {
          ip=indx[i];
          sum=b[ip];
          b[ip]=b[i];
          if (ii)
               for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
          else if (sum) ii=i;
          b[i]=sum;
     }
     for (i=n;i>=1;i--) {
          sum=b[i];
          for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
          b[i]=sum/a[i][i];
     }
  return(0);
}
 


# define  PI         3.141592653589793238462643


//compute determinant of d-dimensional matrix A
double det(double **A)
{
	double det2, det_dmat(); 
	
	int ERROR; 	
	
	if(d==2) 
	{
		det2=(A[1][1]*A[2][2] - A[1][2]*A[1][2]);				
	}
	else if(d==3)
	{
		det2=(A[1][1]*(A[2][2]*A[3][3]   - A[2][3]*A[2][3])  + A[1][2]*(A[2][3]*A[3][1]  -A[2][1]*A[3][3])  + A[1][3]*(A[2][1]*A[3][2]   -  A[2][2]*A[3][1]) );
	}
	else
	{
		det2=det_dmat(A,d, &ERROR);
		//printf("det A=%lf\n", det2);				
	}
			
	if(det2 <=0.0) 
	{	
		printf("#det2 <=0!\n");	
		
		//print  A matrix
		/*printf("#A:=Matrix(%d,%d, [",d,d);
		for(l=1;    l<=d;  l++)
   	{
   			for(s=1;    s<=d;  s++)
   			{
				
					printf("%lf,\t", A[l][s] );
				}
		}		
		printf("]);\n");*/
		

	
	}
	
	return det2;
	
}


void compute_parameters(double **x, int *Mu, int *M, double **m, double **Q, double **P,  double *Det,  double **A,  double **B)
{
	int mu, j, i, l,s;
	
	double det();
	
	
	//compute sizes, mean, covariance, precision  matrix, det. of cov. of clusters	
	
	//prepare arrays 
	for(mu=1;   mu<=K;     mu++)   
    {
		M[mu]=0;
		
		//erase mean	   		
   		for(j=1;    j<=d;  j++) 
   														m[mu][j]=0.0;
   									
  		//erase cov. and precision matrix 
   		for(l=1;    l<=d;  l++)
   		{  
   			Q[(mu-1)*d+l][l]=0.0;
   			P[(mu-1)*d+l][l]=0.0;
   		}
				
		for(l=1;    l < d;  l++)
		{
			for(s=l+1;    s<=d;  s++)
			{
				Q[(mu-1)*d+l][s]=0.0;
				Q[(mu-1)*d+s][l]=0.0;
				P[(mu-1)*d+l][s]=0.0;
				P[(mu-1)*d+s][l]=0.0;							
			}			
		}
   									
   										
   }
	
	//compute covariance matrix
	for(i=1;    i<=N;  i++)
	{
		
		M[Mu[i]]+=1;
		
		for(j=1;    j<=d;  j++)
						m[Mu[i]][j]+=x[i][j];
						
		for(l=1;    l<=d;  l++)  Q[(Mu[i]-1)*d+l][l] +=(x[i][l]*x[i][l]);
				
		for(l=1;    l < d;  l++)
		{
				for(s=l+1;    s<=d;  s++)
																Q[(Mu[i]-1)*d+l][s] +=(x[i][l]*x[i][s]);	
		}
						
	}
	
   for(mu=1;   mu<=K;     mu++)   
   {
			
		for(j=1;    j<=d;  j++) 
   									m[mu][j]/=((double)M[mu]);	
		
		for(l=1;    l<=d;  l++)
		{
			Q[(mu-1)*d+l][l]/=((double)M[mu]);
			Q[(mu-1)*d+l][l]-=(m[mu][l]*m[mu][l]);
		
		}		
				
		for(l=1;    l < d;  l++)
		{
				for(s=l+1;    s<=d;  s++)
				{
					Q[(mu-1)*d+l][s]/=((double)M[mu]);
   					Q[(mu-1)*d+l][s]-=(m[mu][l]*m[mu][s]);
						
				}

		}
										
   }
	
   //compute determinats and inverses of  covariance matrices, i.e. precision matrices
 	for(mu=1;   mu<=K;     mu++) 
	{
  		for(l=1;    l<=d;  l++) A[l][l]=     Q[(mu-1)*d+l][l];
				
		for(l=1;    l < d;  l++)
		{
				for(s=l+1;    s<=d;  s++)
				{
					A[l][s]=Q[(mu-1)*d+l][s];
					A[s][l]=A[l][s];				
				}			
		}
		
		//compute determinant 
		Det[mu]=det(A);

		//compute precision matrix
		
		//matrix B is inverse of A
		inv_dmat(B, A,d);
  
  		for(l=1;    l<=d;  l++)   
  											P[(mu-1)*d+l][l]= B[l][l];
				
		for(l=1;    l < d;  l++)
		{
				for(s=l+1;    s<=d;  s++)
				{
					P[(mu-1)*d+l][s]=B[l][s];
					//P[(mu-1)*d+s][l]=B[l][s];
					
				}
															
		}
  
  } 
 		
}


int max(int i, double **x, double **m, double **Q, double **P,  double *Det)
{
	int mu, l, s, ln_max_mu;
	
	double SUM, ln_max;
	
	mu=1; 
	
	//compute trace
	SUM=0.0;
		
	//off diagonal  contr. 
	for(l=1;    l < d;  l++)
	{
		for(s=l+1;    s<=d;  s++)
                                    SUM+=(P[(mu-1)*d+l][s]*(x[i][l]-m[mu][l])*(x[i][s]-m[mu][s]));
	}
		
	SUM*=2.0;
		
	//diagonal  contr. 
	for(l=1;    l<=d;  l++)  SUM+=(P[(mu-1)*d+l][l]*(x[i][l]-m[mu][l])*(x[i][l]-m[mu][l]));	
	//end compute trace
		
	SUM*=(-0.5);
		
	SUM-=(0.5*log(Det[mu]));
		
	ln_max=SUM;
		
	ln_max_mu=1;
		
		/*if(i==1) 
		{
			printf("mu=%d logP=%lf\n",mu, SUM);
		}*/
			
	for(mu=2;   mu<=K;     mu++)
	{
   		//compute trace
		SUM=0.0;
		
		//off diagonal  contr. 
		for(l=1;    l < d;  l++)
		{
				for(s=l+1;    s<=d;  s++)
																SUM+=(P[(mu-1)*d+l][s]*(x[i][l]-m[mu][l])*(x[i][s]-m[mu][s]));			
		}
		
		SUM*=2.0;
		
		//diagonal  contr. 
		for(l=1;    l<=d;  l++)  SUM+=(P[(mu-1)*d+l][l]*(x[i][l]-m[mu][l])*(x[i][l]-m[mu][l]));	
		//end compute trace
		
		
		SUM*=(-0.5);
		
		SUM-=(0.5*log(Det[mu]));
		
		/*if(i==1) 
		{
			printf("mu=%d logP=%lf\n",mu, SUM);
		}*/
		
		if(SUM>ln_max) 
		{
			ln_max=SUM;
			ln_max_mu=mu;		
		}
			
	
	}
    
	return ln_max_mu;
		
}


double compute_F(double **x, int *Mu,  double **m, double **Q, double **P, double *Det)
{
	int i, l, s;
	
	double SUM, F;
	
	F=0.0;	
	
	for(i=1;   i<=N;     i++)
	{
		
		//compute trace
		/*SUM=0.0;		
		for(l=1;    l<=d;  l++)
   		{
   			for(s=1;    s<=d;  s++)
   								SUM+=(P[(Mu[i]-1)*d+l][s]*(x[i][l]-m[Mu[i]][l])*(x[i][s]-m[Mu[i]][s]));			
   			
		}*/
		
		//new code 
		//compute trace
		SUM=0.0;
		
		for(l=1;    l < d;  l++)
		{
				for(s=l+1;    s<=d;  s++)
																SUM+=(P[(Mu[i]-1)*d+l][s]*(x[i][l]-m[Mu[i]][l])*(x[i][s]-m[Mu[i]][s]));		
		}
		
		SUM*=2.0;
		
		for(l=1;    l<=d;  l++)  SUM+=(P[(Mu[i]-1)*d+l][l]*(x[i][l]-m[Mu[i]][l])*(x[i][l]-m[Mu[i]][l]));	
				
		
		
		
		
		SUM+=log(Det[Mu[i]]);
		
		SUM*=0.5;
		
		F+=SUM; 
		
		
	
	}
	
	return (  (F/((double)N)) +(log(2.0*PI)*0.5*((double)d))    );
	
	
}

void init_clusters(int *empty, int *Mu, int *M)
{
	int i,mu, nu; 
	double drand48();

	if(K>1) 
	{   

		*empty=1;   	
   	
		while(*empty) 
		{   	
			*empty=0;			
			   		
   			for(mu=1; mu<=K; mu++) M[mu]=0; 
		
   			for(i=1;    i<=N;  i++)   
  			{
   	
  					//select nu randomly from 1,...,K
       			nu=(int)(drand48()*K)+1;
       
       			Mu[i]=nu;
       			
       			M[nu]+=1;      		
   			}
   	
   		//check for  empty or small clusters
			for(mu=1; mu<=K; mu++)
			{ 
		
				if (M[mu]<d) *empty=1; 
		
			}	
		}
		//printf("empty=%d\n",empty);
  }
  
  else
  {
   		for(i=1;    i<=N;  i++)  Mu[i]=1;
  }
}

void transform_data( double (*f)(double), double scale, int column, double **x)
{
	int i;
	
	double z; 

	for(i=1;   i<=N;     i++)
	{
		z=(x[i][column]/scale);
		x[i][column]=f(z);	
	}
		
}


int main()
{
	
	FILE *fopen(),*fp, *fp1;

	char name[50],	file_name1[50], file_name2[50];
	
	int i,  j, seed,     t_max, mu, *ivector(), *Mu, *Mu_min, *M,  empty, K1, K2, r, rest, cntr;

	double **dmatrix(), **x, **Q, **P, **m, **A, **B,  *dvector(), *Det, F_min, F=0.0, F1, dF;
	
	//parameter to control a precision of the population dynamics	
	double delta=1.0e-12;    

   //read parametres from the terminal
   scanf("%d %d %d %d %d %d %d  %s", &seed, &N, &d, &K1,  &K2, &rest, &t_max,  name); //input seed for random number generator, sample-size N, dimension of data d, range  for assumed number of clusters K1..K2, number of random restarts rest, maximum number of sweeps over population  t_max, and data-file name. 

   //write parametres to the terminal
   printf("#seed=%d\t N=%d\t  d=%d\t  K1=%d\t  K2=%d\t rest.=%d\t t_max=%d\t  data-file=%s\n", seed, N, d, K1, K2, rest, t_max, name);


  //define matrices and vectors
	x=dmatrix(1,N,1,  d);//matrix to store N (sample size) d-dimensional  points 
	
	Mu=ivector(1,N);// vector to store mu_i=argmax_{\mu} \log Norm(x_i| mean_\mu, Cov_\mu), where x_i is data-point and  mean_\mu, Cov_\mu are mean and covariance of data cluster \mu
	
	Mu_min=ivector(1,N);// vector  to store mu_i corresponding to a local minimum of -log-likelihood F
	
	A=dmatrix(1,d,1, d);//auxiliary matrix
	
	B=dmatrix(1,d,1, d); //auxiliary matrix
	
	m=dmatrix(1,K2,1,  d);//matrix of mean vectors of clusters
	
	Q=dmatrix(1,(K2*d),1,  d);//matrix of covariance matrices of clusters
	
	P=dmatrix(1,(K2*d),1,  d);//matrix of precision  matrices of clusters
	
	Det=dvector(1,K2);// vector of determinants of covariance matrices of clusters
	
	M=ivector(1,K2);// vector of sizes of clusters
	
	//init. data-matrix x
	for(i=1;	i<=N;	i++) {for (j=1; j<=d; j++) x[i][j]=0.0;} 

	load_data(name, x);
	
	//applay log, or any other function, such as asinh, which takes double (divided by scale parameter) and returns double,   transformation to all columns or to selected columns only
	/*for (j=1; j<=d; j++)
							transform_data(log, 1.0, j, x);*/
	
	//print data
	/*printf("data=\n");
	for(i=1;	i<=N;	i++) {for (j=1; j< d; j++) printf("%lf\t",x[i][j]); printf("%lf\n",x[i][d]);  }
	printf("\n");*/
	
	//set random generator
	srand48((long int)seed);  
	
	//loop over cluster sizes K from K1 to K2
	for(K=K1;	K<=K2;	K++)
	{
		//open file to save results of each restart 
		snprintf(file_name1, sizeof file_name1, ".K%d.F.csv",K);
		
		strcpy(file_name2,name);
		
		strcat(file_name2, file_name1);
				
		fp1=fopen(file_name2,"w");
			
		//write column labels
		fprintf(fp1, "time\t"); 
		for(mu=1; mu<=K; mu++) 
																	fprintf(fp1, "M%d\t", mu); 
		fprintf(fp1, "F\n");
					
		//for each K compute F min by random restarts of population dynamics 
		if(K>1)
		{
			//assign N data-points to K clusters randomly and uniformly
			init_clusters(&empty, Mu, M);
			
			//compute parameters of Gaussians
			compute_parameters(x, Mu, M, m, Q,  P,  Det,  A,  B);
			
			//print headers for parameters
			printf("#t\t"); for(mu=1; mu<=K; mu++) printf("M%d\t", mu); printf("F\n");

			printf("%d\t", 0);
			
			//print out cluster sizes
			for(mu=1; mu<=K; mu++) printf("%d\t", M[mu]); 
			
			//compute -log-likelihood 
			F= compute_F(x, Mu,  m, Q, P, Det); //printf("%lf\n",F );
			
			//print out F
			printf("%lf\n",F );
			
			//use F of a random partition to set init. value of  min -log-likelihood F
			F_min=F;
			
			for(r=1;	r<=rest;	r++)
			{
				
 				cntr=1;

				dF=1.0;
				
				
				/*do the gradient descent until the  -log-likelihood F is not changing within the accuracy delta.  Also this loop is terminated
				 when the number of sweeps over the population exceeds t_max or an empty cluster is created */
				while (((dF>delta) && (cntr<=t_max))&& (!empty))
				{
	
					F1=F; //memorize F	
						
					//update population
                    #pragma omp parallel for private(i) // shared variables are kept the same value on all cores, private values can differ
					for( i=1; i<=N; i++)
                                        Mu[i]=max(i, x, m, Q, P, Det);
					
		
					//compute parameters
					compute_parameters(x, Mu, M, m, Q,  P,  Det,  A,  B);
		
					//check for  empty or small clusters
					for(mu=1; mu<=K; mu++)
                                            if (M[mu]<d) empty=1; 
																		
					//print out observables		
					printf("%d\t", cntr); 
		
					for(mu=1; mu<=K; mu++) printf("%d\t", M[mu]); 
																		
					if(!empty) 
					{
						//if there are no clusters smaller than d then compute neg. -log-likelihood F		 	
		 				F = compute_F(x, Mu,  m, Q, P, Det);
		 				
		 				printf("%lf\n", F  );	
					}
		
					//compute energy difference
					dF=fabs(F-F1); 
		
					cntr+=1;
		
				}
	
				if(F<F_min)
				{
					F_min=F;		
			
  					//save  population corresponding to min -log-likelihood F 
  					for(i=1;    i<=N;  i++) 
  																Mu_min[i]=Mu[i];	
  				}
   			
   			
   				//save time, cluster sizes and -log-likelihood
				fprintf(fp1, "%d\t",cntr); 
				for(mu=1; mu<=K; mu++) fprintf(fp1, "%d\t", M[mu]); 
				fprintf(fp1, "%lf\n", F  );
				
				//prepare for next random restart				
				
				//assign N data-points to K clusters randomly and uniformly
				 init_clusters(&empty, Mu, M);
			
				//compute parameters of Gaussians
				compute_parameters(x, Mu, M, m, Q,  P,  Det,  A,  B);
				
				//compute -log-likelihood 
				F= compute_F(x, Mu,  m, Q, P, Det); //printf("%lf\n",F );
				
				//print out observables		
				printf("%d\t", 0); 
		
				for(mu=1; mu<=K; mu++) printf("%d\t", M[mu]); 
				
				printf("%lf\n", F  );
			
			
	
			} //end of rest. loop
		
			fclose(fp1);
		
		}
		
		//if K is not K>1 then assume that K=1 
		else 
		{
			//assign N data-points to K clusters randomly and uniformly
			 init_clusters(&empty, Mu, M);
			
			//compute parameters of Gaussians
			compute_parameters(x, Mu, M, m, Q,  P,  Det,  A,  B);
			
			//compute -log-likelihood  for K=1
			F_min= compute_F(x, Mu,  m, Q, P, Det); 
			
			//save  population corresponding to min -log-likelihood F 
  			for(i=1;    i<=N;  i++) 
  																Mu_min[i]=1;
  			//save time, cluster sizes and -log-likelihood
			fprintf(fp1, "%d\t%d\t%lf\n",0,1,F_min); 
			
			fclose(fp1);
		}
		
		//save data clustering corresponding to min -log-likelihood F  over restarts 
		snprintf(file_name1, sizeof file_name1, ".K%d.clusters.csv",K);
		
		strcpy(file_name2,name);
		
		strcat(file_name2, file_name1);
				
		fp=fopen(file_name2,"w");
	
		for(i=1;    i<=N;  i++)
  		{
					for(j=1;   j<=d;     j++)  fprintf(fp,"%lf\t", x[i][j]); 
					fprintf(fp,"%d\n", Mu_min[i]);
			 
  		} 
 		fclose(fp);
		
}//end of K loop

//free memory
      
   free_dmatrix(x,1,N,1,d);
   
   free_ivector(Mu,1,N);
   
   free_ivector(Mu_min,1,N);
     
   free_dmatrix(A,1,d,1,d);
   
   free_dmatrix(B,1,d,1,d);
   
   //
   free_dmatrix(m,1,K2,1,d);
   
   free_ivector(M,1,K2);
	 
   free_dmatrix(Q,1,(K2*d),1,d);
   
   free_dmatrix(P,1,(K2*d),1,d);
   
   free_dvector(Det,1,K2);
      
	return 0;
}



