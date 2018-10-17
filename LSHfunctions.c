#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "ATD.h"
#include "LSHfunctions.h"

#define W 4

deiktis** createHashTablesHamming(int K, int L)
{
    int i, j, hashTableSize;
    deiktis** ht;
    
    hashTableSize = mypow(2,K); /* giati i pow mas evgaze 2^4 = 15.999 kapoies fores, opote ftiaksame diki mas */
    
    ht = malloc(L*sizeof(deiktis*)); 
    if (ht == NULL)
    {
        printf("Malloc failed...Programme will now be terminated\n");
        exit(-1);
    }
    
    for (i = 0; i < L; i++)
    {
        ht[i] = malloc(hashTableSize * sizeof(deiktis));
        if (ht[i] == NULL)
        {
            printf("Malloc failed...Programme will now be terminated\n)");
            exit(-1);
        }
        
        for (j = 0; j < hashTableSize; j++)
            ht[i][j] = NULL;
    }
    return ht;
}


deiktis** createHashTablesEuklidean(int N, int L)
{
    int i, j, hashTableSize;
    deiktis** ht;
    
    hashTableSize = N;
    
    ht = malloc(L*sizeof(deiktis*));
    if (ht == NULL)
    {
        printf("Malloc failed...Programme will now be terminated\n");
        exit(-1);
    }
    
    for (i = 0; i < L; i++)
    {
        ht[i] = malloc(hashTableSize * sizeof(deiktis));
        if (ht[i] == NULL)
        {
            printf("Malloc failed...Programme will now be terminated\n)");
            exit(-1);
        }
        
        for (j = 0; j < hashTableSize; j++)
            ht[i][j] = NULL;
    }
    return ht;
}

double** createVector(int K, int dimensions)
{
    int i, j, t;
    double **hv, rand_value;
    
    hv = malloc(K*sizeof(double*));
    if (hv == NULL)
    {
        printf("Malloc failed...Programme will now be terminated");
        exit(-1);
    }
    for (i = 0; i < K; i++)
        hv[i] = malloc(dimensions*sizeof(double));
        
        
    for (i = 0; i < K; i++)
    {   
        for (j = 0; j < dimensions; j++)
            kanonikh((&hv[i][j]), &rand_value);
        }
    return hv;
}

double* createTaf(int K)
{
	int i;
	double* T;
	T = malloc(K*sizeof(double));
	if (T == NULL)
    {
        printf("Malloc failed...Programme will now be terminated");
        exit(-1);
    }
    
	for (i = 0; i < K; i++)
	 	T[i] = omoiomorfi(0, W);
	 	
	return T;	
}

int** createR(int K, int L)
{	
	int i, j;
	int** R;
	R = malloc(L*sizeof(int*));
	if (R == NULL)
    {
        printf("Malloc failed...Programme will now be terminated");
        exit(-1);
    }
    for (i = 0; i < L; i++)
        R[i] = malloc(K*sizeof(int));
    
	for (i = 0; i < L; i++)
		for (j = 0; j < K; j++)
	 	R[i][j] = omoiomorfi(0, 1000);
	 	
	return R;	
}

int H(double* p, double** v, double* t, int dimensions, int row)
{
	double sum = 0;
	int i;
	for(i=0; i<dimensions; i++)
		sum = sum + p[i]*v[row][i];
		
	sum = (sum + t[row])/W;
	
	return (int)sum;
}

int Fi(double* p, double** V, double* T, int** R, int dimensions, int K, int TableSize, int which_hashTable)
{
	int i, sum = 0;
    long M = mypow(2, 30) - 5;
    
    for (i = 0; i < K; i++)
    {
        sum = sum + R[which_hashTable][i]*H(p, V, T, dimensions, i);
    }
    if(sum>0)
    	sum = sum%M;
    else
    	sum = M - sum;
	
	return (sum%TableSize);
}


int mypow(int vasi, int ekthetis)
{
	int i, x = 1;
	for(i=0; i<ekthetis; i++)
		x = x*vasi;
	
	return x;
}

int hashThesi(int *Hvectors, int K, double* item)
{
	int i, varos, sum;
	sum =0; 
	varos =  mypow(2, K);
	
	for(i=0; i<K; i++)
	{
		varos = varos/2;
		if(item[Hvectors[i]] == 1.0)
			sum = sum + varos;
	}
	return sum;
}
 
double omoiomorfi(int N, int M)
{
	double d;
	return d = M + (rand() / (RAND_MAX + 1.0))* (N-M);
}

void kanonikh(double* x1, double* x2)
{
	double y1, y2, r;

	do{
		y1 = omoiomorfi(-1, 1);
		y2 = omoiomorfi(-1, 1);
		r = (y1*y1) + (y2*y2);
		
		if(r < 1)
		{
			*x1 = y1*sqrt ((-2*log(r*r))/(r*r));
			*x2 = y2*sqrt ((-2*log(r*r))/(r*r));
			return;
		}
	}while(1);
}

int hammingDistance(double *item1, double *item2, int dimension)
{
	int i, sum = 0;
	for(i=0; i<dimension; i++)
		if(item1[i] != item2[i])
			sum++;	
	
	return sum;
}

double eucledianDistance(double *item1, double *item2, int dimension)
{
	int i;
	double sum = 0;
	for(i=0; i<dimension; i++)
			sum = sum + (item1[i]-item2[i])*(item1[i]-item2[i]);	
	
	return sqrt(sum);
}

double manhattanDistance(double *item1, double *item2, int dimension)
{
	int i;
	double sum = 0;
	for(i=0; i<dimension; i++)
			sum = sum + fabs(item1[i]-item2[i]);	
	
	return sum;
}

double cosineDistance(double *item1, double *item2, int dimension)
{
	int i = 0;
	double ar = 0.0, par1 = 0.0, par2 = 0.0;
	for(i=0; i<dimension; i++)
	{
		ar = ar + item1[i]*item2[i];
		par1 = par1 + item1[i]*item1[i];
		par2 = par2 + item2[i]*item2[i];
	} 
	
	return (ar / (sqrt(par1)*sqrt(par2)) );
}

int DBH_query_thesi(double** array, double* V, int dimensions, int L, int K, int* x1, int *x2, double* t1, int which_L)
{
	//enters_in = DBH_query_thesi(array[QUERY_ID], dimensions, L, K, x1, x2, t, j);
	
	int varos, sum, z, x1_trexon, x2_trexon;
	double hx1x2 = 0.0;
	
	varos =  mypow(2, K); sum = 0;	
	for(z=0; z<K; z++)
	{
    	// gia ton i hashtable, gia to item j, i Z-osti h
        x1_trexon = x1[z];
        x2_trexon = x2[z];
        			
    //	printf("\nx1_trexon = %d, x2_trexon = %d", x1_trexon, x2_trexon);
        			
		hx1x2 = (V[x1_trexon] * V[x1_trexon]) + (V[x2_trexon] * V[x2_trexon]);
		hx1x2 = hx1x2 - (array[x2_trexon][x1_trexon] * array[x2_trexon][x1_trexon]);
		hx1x2 = hx1x2 / (2*array[x2_trexon][x1_trexon]);
		
		
		varos = varos/2;		
		if(hx1x2 >= t1[z])
			sum = sum + varos;
	}
	return sum;			
}


