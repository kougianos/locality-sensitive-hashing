#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "ATD.h"
#include "LSHfunctions.h"
#include "LSH.h"

#define MAX_HAMMING_SIZE 70
#define W 4
#define RADIUS 200

/* orismata: inputFile, K, L, queryFile, outputFile */
/* 
	–d DataHamming.csv –q QueryHamming.csv –k 4 -L 5 -o outputH.txt
	–d DataEuclidean.csv –q QueryEuclidean.csv –k 4 -L 5 -o outputV.txt
	–d DistanceMatrix.csv –q QueryDistanceMatrix.csv –k 4 -L 5 -o outputM.txt
	
	–d DataHamming.csv –q QueryHamming.csv -o output.txt (default K=4 kai L=5)
	
	'i kanena orisma kai ta dinei o xristis apo to pliktrologio
*/


int main(int argc, const char * argv[])
{
	int K, L;
	char inFile[50], qFile[50], outFile[50]; 
	
    if (argc == 1) // an den ta dwsei mesw orismatwn
    {
    	printf("Doste K: "); scanf("%d", &K);
    	printf("\nDoste L: "); scanf("%d", &L);
        printf("\nDoste Input File: "); scanf("%s", inFile);
        printf("\nDoste Query File: "); scanf("%s", qFile);
        printf("\nDoste Output File: "); scanf("%s", outFile);
        LSHAlgorithm(K, L, inFile, qFile, outFile);
    }
    else if(argc==7) // 7 argv simainei oti default ta k kai L
    {
       	LSHAlgorithm(4, 5, argv[2], argv[4], argv[6]);
    }
    else if (argc == 11) // ara orizei kai ta K kai L
    {
    	LSHAlgorithm(atoi(argv[6]), atoi(argv[8]), argv[2], argv[4],argv[10]);
	}
    else
    {
        printf("Lathos plithos orismatwn!");
        exit(-1);
	}
	
    return 0;
}

