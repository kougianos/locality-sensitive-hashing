#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

int mypow(int vasi, int ekthetis);
double omoiomorfi(int N, int M);
void kanonikh(double* x1, double* x2);

/* Euclidean */
deiktis** createHashTablesEuklidean(int K, int L);
double** createVector(int K, int dimensions);
double* createTaf(int K);
int** createR(int K, int L);
int H(double* p, double** v, double* t, int dimensions, int row);
int Fi(double* p, double** V, double* T, int** R, int dimensions, int K, int TableSize, int which_hashTable);
double eucledianDistance(double *item1, double *item2, int dimension);
double manhattanDistance(double *item1, double *item2, int dimension);
double cosineDistance(double *item1, double *item2, int dimension);


/* Hamming */
deiktis** createHashTablesHamming(int K, int L);
int hashThesi(int *Hvectors, int K, double* item);
int hammingDistance(double *item1, double *item2, int dimension);

/* Matrix */
int DBH_query_thesi(double** array, double* V, int dimensions, int L, int K, int* x1, int *x2, double* t1, int which_L);
