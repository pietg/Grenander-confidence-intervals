//
//  main.cpp
//  Banerjee-Wellner tyoe confidence intervals
//
//  Created by Piet on 24/03/20.
//  Copyright (c) 2020 Piet. All rights reserved.
//


#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <fstream>
#include <string.h>
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

int i0,m,n,ngrid;
double a,b,mu,lambda,*f,*f1,*uu,a_hyp,*CI;
double *cumw,*cs,*y,*y2;

#define SQR(x) ((x)*(x))

int     compare(const void *a, const void *b);
void    concavemaj(int n, double cumw[], double cs[],  double y[]);
void    sort(int n, double xx[]);
int     compare(const void *a, const void *b);
double  minmax(double lambda);
double  golden(double a1, double b1, double (*f)(double));
double  LRtest(int n, double y2[]);


// [[Rcpp::export]]

List ComputeIntervals(NumericVector x1, double B1)
{
    int     i,j,k,npoints;
    double  L,L0;
    double  B,*grid,*p,*jumploc,step;
    double  *lowbound,*upbound;
    
    B=(double)B1;
    
    n = (int)(x1.size());
    a=0;
    
    printf("\n\n");
        
    ngrid=100;
    npoints=ngrid-10;
    
    uu= new double[n+2];
    y= new double[n+2];
    y2= new double[n+2];

    cumw= new double[n+2];
    cs= new double[n+2];
    f= new double[n+2];
    
    p= new double[n+2];
    jumploc= new double[n+2];
    
    grid= new double[ngrid+2];
    CI= new double[ngrid+1];
    
    upbound  = new double[npoints+1];
    lowbound = new double[npoints+1];
    
    step=B/ngrid;
    for (i=0;i<=npoints;i++)
        grid[i]=i*step;
        
    uu[0]=0;
    uu[n+1]=B;
    
    for (i=0;i<n;i++)
      uu[i+1]=(double)x1[i];
           
    sort(n,uu);
    
    cumw[0]=cs[0]=0;
    
    cumw[n+1]=B;
    cs[n+1]=1;
    
    for (i=1;i<=n;i++)
    {
        cumw[i]=uu[i];
        cs[i]=cs[i-1]+1.0/n;
    }
    
    concavemaj(n+1,cumw,cs,y);
    
    L=0;
    for (i=1;i<=n;i++)
        L += log(y[i]);
    
    j=0;
    
    jumploc[1]=uu[1];
    f[1]=y[1];
    
    for (i=2;i<=n+1;i++)
    {
        if (y[i]<y[i-1])
        {
            j++;
            jumploc[j]=uu[i-1];
            p[j]=y[i-1]-y[i];
            f[j]=y[i-1];
        }
    }
    
    m=j;

    jumploc[m+1]=B;
    f[m+1]=0;
    
    for (k=1;k<=npoints;k++)
    {
        i=1;
        while (uu[i]<=grid[k] && i<n)
            i++;
        
        i0=i;
    
        a=y[i0];
        
        CI[k]=a;
        
        a_hyp=a;
        b=0;
        
        j=0;
        
        lambda=1;
        
        while (b<2.26916)
        {
            j++;
            a_hyp=a+j*0.001;
            
            lambda = golden(1.0,10,minmax);
            
            L0 = LRtest(n,y2);
        
            b=2*(L-L0);
        }
        
        upbound[k]=a_hyp;

        j=0;
        
        a=y[i0];
        
        a_hyp=a;
        
        lambda=1;
        
        b=0;
        
        if (a-0.001<=0)
            lowbound[k]=0;
        else
        {
            while (b<2.26916 && a-j*0.001>0.001)
            {
                j++;
                a_hyp=a-j*0.001;
                
                lambda = golden(0.1,1.0,minmax);
                
                L0 = LRtest(n,y2);
                
                b=2*(L-L0);
                
            }
            lowbound[k]=fmax(0,a_hyp);
        }
        Rcout  << setw(10) << k << setprecision(6) <<  setw(15) << lowbound[k] << setprecision(6) <<  setw(15) << CI[k] << setprecision(6) <<  setw(15) << upbound[k]<< setprecision(6) <<  setw(15)<< std::endl;
    }

    NumericMatrix out1 = NumericMatrix(m,2);
    
    for (i=0;i<m;i++)
    {
      out1(i,0)=jumploc[i+1];
      out1(i,1) = f[i+1];
    }
    
    NumericMatrix out2 = NumericMatrix(npoints,4);
    
    for (i=0;i<npoints;i++)
    {
        out2(i,0)=grid[i+1];
        out2(i,1)=lowbound[i+1];
        out2(i,2)=CI[i+1];
        out2(i,3)=upbound[i+1];
    }

    Rcout << "Making output list" << std::endl;
    
    // make the list for the output, containing the MLE and the confidence intervals
    
    List out = List::create(Rcpp::Named("MLE")=out1,Rcpp::Named("CI")=out2);
    
    Rcout << "Freeing memory" << std::endl;
    
    // free memory
    
    delete[] uu; delete[] p, delete[] jumploc; delete[] grid; delete[] CI;
    delete[] f; delete[] cumw; delete[] cs;
    delete[] y; delete[] y2; delete[] lowbound; delete[] upbound;
    
    return out;
    
}

int compare(const void *a, const void *b)
{
    double x = *(double*)a;
    double y = *(double*)b;
    
    if (x < y)
        return -1;
    if (x > y)
        return 1;
    return 0;
}


void sort(int n, double xx[])
{
    int i;
    double *obs;
    
    obs = new double[n];
    
    for (i=0;i<n;i++)
        obs[i]=xx[i+1];
    
    qsort(obs,n,sizeof(double),compare);
    
    for (i=1;i<=n;i++)
        xx[i]=obs[i-1];
    
    delete[] obs;
}

double    LRtest(int n, double y2[])
{
    int     i;
    double  sum;
    
    for (i=1;i<=n;i++)
        cumw[i] = lambda*uu[i];
    
    for (i=1;i<=n;i++)
    {
        cs[i] = cs[i-1]+1.0/n;
        if (i==i0)
            cs[i] += lambda-1;
    }
    
    cs[n+1]=cs[n];
    
    concavemaj(n+1,cumw,cs,y2);
    
    sum=0;
    for (i=1;i<=n;i++)
        sum += log(y2[i]);
    
    return sum;
}

double minmax(double lambda)
{
    int i,j,k;
    double sum,min,max;
    
    min=1000;
    
    for (i=1;i<=i0;i++)
    {
        max=0;
        for (j=i0;j<=n;j++)
        {
            sum=0;
            for (k=i;k<=j;k++)
                sum += 1.0/n;
            sum += lambda-1;
            
            sum /= ((uu[j]-uu[i-1])*lambda);
           
            
            if (sum>max)
                max=sum;
        }
        if (max<min)
            min=max;
    }
    return fabs(a_hyp-min);
}


double golden(double a1, double b1, double (*f)(double))
{
    double a,b,eps=1.0e-10;
    
    a=a1;
    b=b1;
    
    double k = (sqrt(5.0) - 1.0) / 2;
    double xL = b - k*(b - a);
    double xR = a + k*(b - a);
    
    while (b-a>eps)
    {
        if ((*f)(xL)<(*f)(xR))
        {
          b = xR;
          xR = xL;
          xL = b - k*(b - a);
        }
        else
        {
          a = xL;
          xL = xR;
          xR = a + k * (b - a);
        }
    }
    
    return (a+b)/2;
}


void concavemaj(int n, double cumw[], double cs[],  double y[])
{
    int    i,j,m;
    
    y[1] = cs[1]/cumw[1];
    for (i=2;i<= n;i++)
    {
        if (cumw[i]-cumw[i-1]>0)
            y[i] = (cs[i]-cs[i - 1])/(cumw[i]-cumw[i-1]);
        else
            y[i]=0;
        if (y[i-1]<y[i])
        {
            j = i;
            while ((y[j-1] < y[i]) && (j>1))
            {
                j=j-1;
                y[i] = (cs[i]-cs[j-1])/(cumw[i]-cumw[j-1]);
                for (m=j;m<i;m++)    y[m] = y[i];
            }
        }
    }
    
}
