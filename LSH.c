#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
//#include <sys/time.h>
#include "ATD.h"
#include "LSHfunctions.h"

#define MAX_HAMMING_SIZE 70
#define W 4
#define MAXDISTANCE 999999999

void sort(double* array, int arrayLength);

void LSHAlgorithm(int K, int L, char* inFile, char* qFile, char* outFile)
{
	FILE* infile, *outfile, *queryfile;
	struct timeval time1, time2;
    char dummy[30], dataType[30], distance_metric, xar, found, ham[MAX_HAMMING_SIZE], **itemNames;
    char queryItemName[30], queryH[MAX_HAMMING_SIZE];
    int dimensions, j, plithos_items, i,  vec, m, loops, enters_in, **ER, **Hvectors;
	int QUERY_ID, **x1, **x2, z, x1_trexon, x2_trexon, varos, sum, NNid;
    double distance_matrix, *ET, **EV, **array, ***hx1x2, **t1, *h_hold, RADIUS, *queryHtoDouble, minDist, current_dist, sumtime1=0, sumtime2=0, elapsedTime;
	deiktis** HT, temp;
	char* alreadyPrinted;
	
	srand(time(NULL));

    infile = fopen(inFile, "r");
    if(infile == NULL) {printf("Den mporei na anoiksei to arxeio eisodou  -> EXIT"); exit(-1);}
    queryfile = fopen(qFile, "r");
    if(queryfile == NULL) {printf("Den mporei na anoiksei to arxeio query  -> EXIT"); exit(-1);}
    outfile = fopen(outFile, "w");
    if(outfile == NULL) {printf("Den mporei na anoiksei to arxeio output  -> EXIT"); exit(-1);}
    
 	fscanf(infile,"%s %s", dummy, dataType);
 	printf("\n %s  %s", dummy, dataType);

	if(dataType[0] == 'v' || dataType[0] == 'e') // vector
 	{
 		fscanf(infile,"%s %s", dummy, dataType);
 		printf("\n %s  %s", dummy, dataType);
 		distance_metric = dataType[0]; // 'e' or 'm' or 'c'
 		// vriskw diastaseis pinaka (plithos items kai dimensions)
 		
 		 fscanf(infile, " %c", &xar); // me keno, gia na min parei to '\n'
		 for(dimensions = 0; xar!='\n';)
		 {
		 	if(xar == '\t') dimensions++; //oi diastaseis xorizontai metaksi tous me tabs, opote auksanoume dimensions gia kathe /t pou sinantame
		 	fscanf(infile, "%c", &xar);
		 }
 		printf("\n DIMENSIONS = %d", dimensions);
 		
 		
 		for(plithos_items = 0; !feof(infile); plithos_items++)
 		{
 			for(j=0; j<dimensions+1; j++) // onoma item + diastaseis tou
 				fscanf(infile, "%s", dummy);
		}
 		printf("\n plithos_items = %d", plithos_items);
 		
	
 		array = (double**)malloc(plithos_items*sizeof(double*)); //pianoume xoro gia disdiastato pinaka pou tha periexei tis times twn items
 		for(i=0; i<plithos_items; i++)
		 	array[i] = (double*)malloc(dimensions*sizeof(double));
		
		itemNames = (char**)malloc(plithos_items*sizeof(char*)); //pianoume xoro gia pinaka pou tha periexei ta onomata twn items
 		for(i=0; i<plithos_items; i++)
		 	itemNames[i] = (char*)malloc(20*sizeof(char));		
 		
 		fseek(infile, 0, SEEK_SET); // pame to deikti tou arxeiou apostasi 0 apo tin arxi
 		fscanf(infile,"%s %s %s %s", dummy, dummy, dummy, dummy); //diavazoume tis protes 4 sumvoloseires pou den mas endiaferoun
 		for(i=0; i<plithos_items; i++)
 		{
 			fscanf(infile, "%s", itemNames[i]);	//diavazoume to onoma tou item kai to vazoume stin i thesi tou pinaka itemNames
 			for(j=0; j<dimensions; j++)
 				fscanf(infile, "%lf", &array[i][j]); //diavazoume tin diastasi tou ekastote i item kai tin vazoume stin j thesi tou pinaka array
		}			
	//	printf("\n %lf %lf %s", array[1][4], array[999][97], itemNames[65]);
 	}
	else if(dataType[0] == 'h') // hamming
	{
		distance_metric = 'h';
		for(plithos_items = -1; !feof(infile); plithos_items++) //-1 gia to prwto enter, sto telos tis 1is grammis
			fscanf(infile,"%s %s", dummy, ham);
		
		dimensions = strlen(ham); //i deuteri simvoloseira pou diavazoume se kathe grammi einai oi diastaseis tou item
			
		printf("\n plithos_items = %d, dimensions = %d", plithos_items, dimensions);
	
		array = (double**)malloc(plithos_items*sizeof(double*)); //pianoume xoro gia disdiastato pinaka pou tha periexei tis times twn items
 		for(i=0; i<plithos_items; i++)
		 	array[i] = (double*)malloc(dimensions*sizeof(double));
		
		itemNames = (char**)malloc(plithos_items*sizeof(char*)); //pianoume xoro gia pinaka pou tha periexei ta onomata twn items
 		for(i=0; i<plithos_items; i++)
		 	itemNames[i] = (char*)malloc(20*sizeof(char));
		 	
		 	
		fseek(infile, 0, SEEK_SET); // pame to deikti tou arxeiou apostasi 0 apo tin arxi
 		fscanf(infile,"%s %s", dummy, dummy);
 		for(i=0; i<plithos_items; i++)
 		{
 			fscanf(infile, "%s %s", itemNames[i], ham);	//se kathe grammi i proti sumvoloseira einai to onoma tou item kai i deuteri i timi tou
 			for(j=0; j<dimensions; j++)
 				array[i][j] = ham[j] - '0';
		}	
		/*printf("\n plithos_items = %d, dimensions = %d\n", plithos_items, dimensions);	
		for(i=0; i<dimensions; i++)
			printf("%1.0lf", array[3][i]);	*/
	}
	else if(dataType[0] == 'm') // matrix
	{
		distance_metric = 't'; // maTrix
		plithos_items = 0;
		
		for(i=0; !feof(infile); i++)
		{
			fscanf(infile, " %c", &xar);
			if(xar == ',') plithos_items++; //kathe fora pou vriskoume komma auksanoume to plithos items kata 1
			if(xar=='\n') break;
		}
		plithos_items++; //giati to teleutaio item den exei komma
		
		printf("\nDM plithos_items = %d\n", plithos_items);
		array = (double**)malloc(plithos_items*sizeof(double*)); //pianoume xoro gia  disdiastato pinaka pou tha periexei tis apostaseis twn items
 		for(i=0; i<plithos_items; i++) array[i] = (double*)malloc(plithos_items*sizeof(double));
		itemNames = (char**)malloc(plithos_items*sizeof(char*)); //pianoume xoro gia pinaka pou tha periexei ta onomata twn items
 		for(i=0; i<plithos_items; i++)itemNames[i] = (char*)malloc(20*sizeof(char));
 		
 		fseek(infile, 0, SEEK_SET); // pame to deikti tou arxeiou apostasi 0 apo tin arxi
		fscanf(infile,"%s %s %s", dummy, dummy, dummy); //oi treis protes sumvoloseires den mas endiaferoun
		
		for(i=0, j=0; i<plithos_items-1;)
		{
			fscanf(infile, " %c", &xar);
			if(xar==',') //an diavasoume komma, simainei oti exoume diavasei ena olokliro itenName
			{
				dummy[j]='\0'; //opote vazoume ton teleutaio xaraktira tou itemName na einai \0 gia na exoume simvoloseira
				strcpy(itemNames[i], dummy); //antigrafi tou itemName pou exoume diavasei mexri tora stin i thesi tou pinaka itemNames
				j = 0; //vazoume to j=0 giati tha ksekinisoume na diavazoume neo onoma
				i++; //pame stin epomeni thesi tou pinaka itemNames
			}
			else //an den diavasoume komma, simainei oti diavazoume xaraktira pou anikei sto onoma tou item
			{
				dummy[j] = xar;
				j++;
			}
		}
		fscanf(infile,"%s", itemNames[plithos_items-1]); //sto teleutaio item den sunantame komma, opote to kanoume me to xeri	
		printf("\n LAST ITEM= %s\n", itemNames[plithos_items-1]);
		
		for(i=0; i<plithos_items; i++)
 		{
 			for(j=0; j<plithos_items; j++)
 				fscanf(infile, "%lf", &array[i][j]);	//gemizoume ton disdiastato pinaka me tis apostaseis ton items
		}	
		//printf("\n %lf %lf %s", array[1][4], array[999][997], itemNames[65]);	
	}
	else // wrong fileType
	{
		printf("\nWrong file type\n Exiting program :-( \n");
		exit(-1);
	}
	
	
	/* -------------------------------  Hamming  ------------------------------- */
	
	if(distance_metric == 'h')
	{	
		Hvectors = (int**)malloc(L*sizeof(int*));
 		for(i=0; i<L; i++) Hvectors[i] = (int*)malloc(K*sizeof(int));
		
		printf("\n K = %d, L = %d, dimensions = %d", K, L, dimensions);
		
		for(vec=0; vec<L;)
		{
	    	for (i = 0; i < K;)		/* gia ta ypoloipa K-1 */
		    {
		        Hvectors[vec][i] = rand() % dimensions; /* dialegw arxika 1 stin tixi */
		        
		        found = 0;
		        for (j = 0; j < i; j++)    /* elegxw na min yparxei idi stis epiloges mou, gia na min parw 2 fores to idio */
		        {
		            if (Hvectors[vec][i] == Hvectors[vec][j]) /* an yparxei idi */
		                found = 1;	/* den proxwraw parakatw kai anebenw sto loop gia tin idia thesi */
		        }
		        if (found == 0)      /* alliws paw stin epomeni thesi */ 
		            i++;
		    }
		    
		    /* prepei episis na ginei elegxos oti ta L dianysmata den einai idia metaksy tous */
		    
		    for (i = 0; i < vec; i++) /* gia kathe ena apo ta proigoumena tou */
		    {
		    	found = 0;
			    for (j = 0; j < K; j++)    /* elegxw na min yparxei idi stis epiloges mou, gia na min parw 2 fores to idio */
			    {
			    	for(m = 0; m < K; m++)
			    	{
			    		if (Hvectors[vec][j] == Hvectors[i][m]) /* an yparxei idi */
			            	found++;	
					}
				}
				if(found == K)
			    	break;
			}
				
			if (found < K)      /* alliws paw stin epomeni thesi */ 
		       vec++;
		}
		
	//	for(vec=0; vec<L; vec++)
	//	{
	//		printf("\n");
	//		for (i = 0; i < K; i++)
	//			printf("%d ", Hvectors[vec][i]);
	//	}	

		/* telos dimiourgias hashvectors */
	
		HT = createHashTablesHamming(K, L); //dimiourgia kenou pinaka hashTables
		
		for(i=0; i<plithos_items; i++)
		{
			for(j=0; j<L; j++)                                         
			{
				enters_in =  hashThesi(Hvectors[j], K, array[i]);
				eisagogi_arxi( &(HT[j][enters_in]), i);
			}
		}

		///////////////////////// QUERY /////////////////////////////////
		fscanf(queryfile, "%s", dummy); //word "Radius:"
		fscanf(queryfile, "%lf", &RADIUS);
		
		queryHtoDouble = (double*)malloc(dimensions*sizeof(double)); //gia na metatrepsei to %s se pinaka %lf
		alreadyPrinted = (char*)malloc(plithos_items*sizeof(char)); // gia na min ektipwsei 2 fores to kathe item
		
		while(!feof(queryfile))
		{
			fscanf(queryfile, "%s", queryItemName);
			if(feof(queryfile))
				break;
			
			fprintf(outfile,"Query: %s", queryItemName);
			fscanf(queryfile, "%s", queryH); // diavazei san string tis coordinates
			
			for(j=0; j<dimensions; j++)
 				queryHtoDouble[j] = queryH[j] - '0'; // tis metatrepei se pinaka double
			
			for(j=0; j<plithos_items; j++) // arxika den exei ektypwthei kanena item san neighbor
				alreadyPrinted[j] = '0';
			
			minDist = MAXDISTANCE; NNid = -1;
			gettimeofday(&time1,NULL);
			for(j=0; j<L; j++)
			{
				enters_in = hashThesi(Hvectors[j], K, queryHtoDouble); // se poia thesi tou hashtable psaxnei
				
				for(temp = HT[j][enters_in]; temp!=NULL; temp = temp->next) // oli ti lista
				{					
					// apostasi apo to kathe item pou vriskei sti lista an einai entos R
					if(hammingDistance(queryHtoDouble, array[temp->indexNo], dimensions) < RADIUS)
					{
						//an exei idi ektypwthei (apo proigoumeno pinaka), proxorame
						if(alreadyPrinted[temp->indexNo] == '1')
							continue;
							
						alreadyPrinted[temp->indexNo] = '1';
						// alliws print sto outfile
						fprintf(outfile,"\n %s", itemNames[temp->indexNo]);
						
						// kai sygkrisi gia min distance
						if(minDist > hammingDistance(queryHtoDouble, array[temp->indexNo], dimensions))
						{
							minDist = hammingDistance(queryHtoDouble, array[temp->indexNo], dimensions);
							NNid = temp->indexNo;
						}
					} 
				}
			}
			gettimeofday(&time2,NULL); 
			if(NNid>-1)// an vrike neighbor(s)
			{
				fprintf(outfile,"\nNearest neighbor: %s", itemNames[NNid]);
				fprintf(outfile,"\ndistanceLSH: %lf", minDist);
			}
			else
			{
				fprintf(outfile,"\nNearest neighbor: -");
			}
			gettimeofday(&time1,NULL);
			//ypologizei (eksantlitika - O(N)) tin pragmatiki min distance
			minDist = MAXDISTANCE; NNid = -1;
			for(j=0; j<plithos_items; j++)
				if(minDist > hammingDistance(queryHtoDouble, array[j], dimensions))
				{
					minDist = hammingDistance(queryHtoDouble, array[j], dimensions);
					NNid = j;
				}
			gettimeofday(&time2,NULL); 
			elapsedTime = (time2.tv_sec - time1.tv_sec) * 1000.0;      /* sec to ms */
      			elapsedTime += (time2.tv_usec - time1.tv_usec) / 1000.0;   /* us to ms */
			//fprintf(outfile, "tLSH: %g", elapsedTime);
			fprintf(outfile,"\nNearest neighbor True: %s", itemNames[NNid]);
			fprintf(outfile,"\ndistanceTrue: %lf\n\n", minDist);
		}
		for(i=0;i<L;i++){free(Hvectors[i]); Hvectors[i]=NULL;} free(Hvectors); //free HVectors
	}
	
	/* -------------------------------  Euclidean  ------------------------------- */
	else if(distance_metric == 'm' || distance_metric == 'e' || distance_metric == 'c')
	{
		HT = createHashTablesEuklidean(plithos_items/4, L);
		EV = createVector(K, dimensions);
		ET =  createTaf(K);
		ER = createR(K, L);
		
		for(i=0; i<plithos_items; i++)
		{
			for(j=0; j<L; j++)
			{
				enters_in = Fi(array[i], EV, ET, ER, dimensions, K, plithos_items/4, j);
				//printf(" enters_in = %d", enters_in);
				eisagogi_arxi( &(HT[j][enters_in]), i);
			}
		}
	
	    /////////// QUERY //////////////
	    fscanf(queryfile, "%s", dummy); //word "Radius:"
		fscanf(queryfile, "%lf", &RADIUS);
		
		queryHtoDouble = (double*)malloc(dimensions*sizeof(double));
		alreadyPrinted = (char*)malloc(plithos_items*sizeof(char));
		
		while(!feof(queryfile))
		{
			fscanf(queryfile, "%s", queryItemName);
			if(feof(queryfile))
				break;
			
			fprintf(outfile,"Query: %s", queryItemName);
			for(j=0; j<dimensions; j++)
 				fscanf(queryfile, "%lf", &(queryHtoDouble[j]));
			
			for(j=0; j<plithos_items; j++)// arxika den exei ektypwthei kanena item san neighbor
				alreadyPrinted[j] = '0';
			
			minDist = MAXDISTANCE; NNid = -1;
			gettimeofday(&time1,NULL);
			for(j=0; j<L; j++) // gia kathe hashtable
 			{		
 				//se poia thesi tou (lista) psaxnei
				enters_in = Fi(queryHtoDouble, EV, ET, ER, dimensions, K, plithos_items/4, j);
				
				// gia kathe stoixeio tis listas
				for(temp = HT[j][enters_in]; temp!=NULL; temp = temp->next)
				{					
					//analogw tin metriki apostasis
					if(distance_metric == 'e')
						current_dist = eucledianDistance(queryHtoDouble, array[temp->indexNo], dimensions);
					else if(distance_metric == 'c')
						current_dist = cosineDistance(queryHtoDouble, array[temp->indexNo], dimensions);
					else// manhattan
						current_dist = manhattanDistance(queryHtoDouble, array[temp->indexNo], dimensions);
				
				//printf("\n ITEM %s distance %lf", itemNames[temp->indexNo], eucledianDistance(queryHtoDouble, array[temp->indexNo], dimensions));
					// an einai entos R
					if(current_dist < RADIUS)
					{	
						// an exei idi ektypwthei (apo proigoumeno hashtable), den to theloume pali
						if(alreadyPrinted[temp->indexNo] == '1')
							continue;
							
						alreadyPrinted[temp->indexNo] = '1';
					
						// alliws ektypwsi
						fprintf(outfile,"\n %s", itemNames[temp->indexNo]);
						// kai elegxos kai enimerwsi gia min distance
						if(minDist > current_dist)
						{
							minDist = current_dist;
							NNid = temp->indexNo;
						}
					} 
				}
			}
			gettimeofday(&time2,NULL); 	
			elapsedTime = (time2.tv_sec - time1.tv_sec) * 1000.0;      /* sec to ms */
      			elapsedTime += (time2.tv_usec - time1.tv_usec) / 1000.0;   /* us to ms */
			//fprintf(outfile, "tLSH: %g", elapsedTime);

			// an vrike estw kai 1 neighbor ektypwnei
			if(NNid>-1)
			{
				fprintf(outfile,"\nNearest neighbor: %s", itemNames[NNid]);
				fprintf(outfile,"\ndistanceLSH: %lf", minDist);
			}
			else
			{
				fprintf(outfile,"\nNearest neighbor: -");
			}
			gettimeofday(&time1,NULL);
			//ypologizei (eksantlitika - O(N)) tin pragmatiki min distance
			minDist = MAXDISTANCE; NNid = -1;
			for(j=0; j<plithos_items; j++)
			{
				if(distance_metric == 'e')
					current_dist = eucledianDistance(queryHtoDouble, array[j], dimensions);
				else if(distance_metric == 'c')
					current_dist = cosineDistance(queryHtoDouble, array[j], dimensions);
				else// manhattan
					current_dist = manhattanDistance(queryHtoDouble, array[j], dimensions);
			
			
				if(minDist > current_dist)
				{
					minDist = eucledianDistance(queryHtoDouble, array[j], dimensions);
					NNid = j;
				}
			}			
			gettimeofday(&time2,NULL); 
			elapsedTime = (time2.tv_sec - time1.tv_sec) * 1000.0;      /* sec to ms */
        		elapsedTime += (time2.tv_usec - time1.tv_usec) / 1000.0;   /* us to ms */
        		//fprintf(outfile, "tLSH: %g", elapsedTime);
			fprintf(outfile,"\nNearest neighbor True: %s", itemNames[NNid]);
			fprintf(outfile,"\ndistanceTrue: %lf\n\n", minDist);
		}
	
		for(i=0; i<K;i++){free(EV[i]);EV[i]=NULL;} free(EV);
		free(ET);
		for(i=0;i<K;i++){free(ER[i]);ER[i]=NULL;} free(ER);
		
	}
	else // if(distance_metric == 't') DBH
	{
		HT = createHashTablesEuklidean(mypow(2,K), L);
		dimensions = plithos_items;
		
		/* dimourgia twn x1 kai x2 gia kathe L */
		x1 = (int**)malloc(L*sizeof(int*));
		x2 = (int**)malloc(L*sizeof(int*));
		
		for(i=0; i<L; i++)
		{
			x1[i] = (int*)malloc(K*sizeof(int));
			x2[i] = (int*)malloc(K*sizeof(int));
		}
		
		printf("\n dimensions = %d", dimensions);
		for(i=0; i<L; i++)	 /* Gia kathenan apo tous L hashtables, dialegw K pairs x1 kai x2 */
		{
		    for (j=0; j<K; j++)
		    {
		        for(;;)
		        {
		        	x1[i][j] = rand()%dimensions;		/* Holds index of item */
		            x2[i][j] = rand()%dimensions;	/* Holds index of item */
		            if(x1[i][j] != x2[i][j])
		            	break;
		        }
		    }
		}
		
		/***** ypologismos h gia kathe L ********/
		hx1x2 = (double***)malloc(L*sizeof(double**));
	 	/* For L hash functions g */
        for (i=0; i<L; i++)
        {
			hx1x2[i] = (double**)malloc(dimensions*sizeof(double*));
			for (j=0; j<dimensions; j++)
		 	{
		 		hx1x2[i][j] = (double*)malloc(sizeof(double)*K);
    
	        	for (z=0; z<K; z++)
    			{
    			//	printf("\n %d-%d-%d", i, j,z);
    				// gia ton i hashtable, gia to item j, i Z-osti h
        			x1_trexon = x1[i][z];
        			x2_trexon = x2[i][z];
        			
        		//	printf("\nx1_trexon = %d, x2_trexon = %d", x1_trexon, x2_trexon);
        			
					hx1x2[i][j][z] = (array[j][x1_trexon] * array[j][x1_trexon]) + (array[j][x2_trexon] * array[j][x2_trexon]);
					hx1x2[i][j][z] = hx1x2[i][j][z] - (array[x2_trexon][x1_trexon] * array[x2_trexon][x1_trexon]);
					hx1x2[i][j][z] = hx1x2[i][j][z] / (2*array[x2_trexon][x1_trexon]);
        		}
           	}
		}
        
        /************** ypologismos tou t1 (to t2 = +oo ) *************/
        t1 = (double**)malloc(L*sizeof(double*));
		for(i=0; i<L; i++)    
            t1[i] = (double*)malloc(K*sizeof(double));
            
    	h_hold = (double*)malloc(dimensions*sizeof(double));
            
	    for (i = 0; i < L; i++)
	    {
	        for (j=0; j<K; j++)
	        {
	       
	        	for(z=0; z<dimensions; z++)
	        	{
					
	        		h_hold[z] = hx1x2[i][z][j];
	        	}
	       
	        	sort(h_hold, dimensions);
	        	
	            t1[i][j] = h_hold[dimensions/2];
	        }
	    }

    	for (i=0; i<L; i++) //gia kathe pinaka
        {
			for (j=0; j<dimensions; j++) //kai kathe item
		 	{
    			//ypologizw se poia thesi tou i-ostou hashtable prepei na paei
    			
    			varos =  mypow(2, K); sum = 0;	
    			for(z=0; z<K; z++)
				{
					varos = varos/2;
					
					if(hx1x2[i][j][z] >= t1[i][z])
						sum = sum + varos;
				}			
    			
    			enters_in = sum;
			//	printf("\n (%d,%d)  ---->> enters_in = %d", i, j, enters_in);
				eisagogi_arxi( &(HT[i][enters_in]), j);
    		}
    	}
    	
    	////////////////////      QUERY      ///////////////////
    
    
    	fscanf(queryfile, "%s", dummy); //word "Radius:"
		fscanf(queryfile, "%lf", &RADIUS);
		
		queryHtoDouble = (double*)malloc(dimensions*sizeof(double));
		alreadyPrinted = (char*)malloc(plithos_items*sizeof(char));
		
		while(!feof(queryfile))
		{
			fscanf(queryfile, "%s", queryItemName);
			if(feof(queryfile))
				break;
			
			fprintf(outfile,"Query: %s", queryItemName);
			for(j=0; j<dimensions; j++)
 				fscanf(queryfile, "%lf", &(queryHtoDouble[j]));
			
			// arxika den exei kanena item ws neihgbor
			for(j=0; j<plithos_items; j++)
				alreadyPrinted[j] = '0';
			gettimeofday(&time1,NULL);
			minDist = MAXDISTANCE; NNid = -1;
			for(j=0; j<L; j++) // gia kathe hashtable
			{		
				//vriskei se poia thesi (lista) na psaxei
				enters_in = DBH_query_thesi(array, queryHtoDouble, dimensions, L, K, x1[j], x2[j], t1[j], j);
				
				// akoluthei oli ti lista
				for(temp = HT[j][enters_in]; temp!=NULL; temp = temp->next)
				{					
					current_dist = queryHtoDouble[temp->indexNo];
				
					// an to item einai entos R
					if(current_dist < RADIUS)
					{	
						// kai den exei idi ektypwthei
						if(alreadyPrinted[temp->indexNo] == '1')
							continue;
							
						alreadyPrinted[temp->indexNo] = '1';
					
						//to ektipwnei
						fprintf(outfile,"\n %s", itemNames[temp->indexNo]);
						//elegxos kai enimerwsi gia NN
						if(minDist > current_dist)
						{
							minDist = current_dist;
							NNid = temp->indexNo;
						}
					} 
				}
			}
			gettimeofday(&time2,NULL); 
			elapsedTime = (time2.tv_sec - time1.tv_sec) * 1000.0;      /* sec to ms */
      			elapsedTime += (time2.tv_usec - time1.tv_usec) / 1000.0;   /* us to ms */
			//fprintf(outfile, "tLSH: %g", elapsedTime);

			// an vrike neighbor(s) ektypwnei 
			if(NNid>-1)
			{
				fprintf(outfile,"\nNearest neighbor: %s", itemNames[NNid]);
				fprintf(outfile,"\ndistanceLSH: %lf", minDist);
			}
			else
			{
				fprintf(outfile,"\nNearest neighbor: -");
			}
			gettimeofday(&time1,NULL);
			//ypologizei (eksantlitika - O(N)) tin pragmatiki min distance
			minDist = MAXDISTANCE; NNid = -1;
			for(j=0; j<plithos_items; j++)
				if(minDist > queryHtoDouble[j])
				{
					minDist = queryHtoDouble[j];
					NNid = j;
				}
			gettimeofday(&time2,NULL); 
			elapsedTime = (time2.tv_sec - time1.tv_sec) * 1000.0;      /* sec to ms */
        		elapsedTime += (time2.tv_usec - time1.tv_usec) / 1000.0;   /* us to ms */
        		//fprintf(outfile, "tLSH: %g", elapsedTime);
			fprintf(outfile,"\nNearest neighbor True: %s", itemNames[NNid]);
			fprintf(outfile,"\ndistanceTrue: %lf\n\n", minDist);
		}
    
    	for(i=0;i<L;i++) //free x1,x2
    	{
    		free(x1[i]); x1[i]=NULL;
    		free(x2[i]); x2[i]=NULL;
    	}
    	free(x1);free(x2);
    	for(i=0;i<L;i++) //free hx1x2
    	{
    		for(j=0;j<dimensions;j++)
    		{
    			free(hx1x2[i][j]); hx1x2[i][j]=NULL;
    		}
    		free(hx1x2[i]); hx1x2[i]=NULL;
    	}
    	free(hx1x2);
    	free(h_hold);
    	for(i=0;i<L;i++){free(t1[i]); t1[i]=NULL;} free(t1);  //free t1
    	
	}
	for(i=0; i<plithos_items; i++)                //free itemNames, array
	{
		free(itemNames[i]); free(array[i]);
		itemNames[i] = NULL; array[i] = NULL;
	}
	free(itemNames); free(array);
	itemNames = NULL; array = NULL;	
}


void sort(double* array, int arrayLength)
{
    int i, j;
    double temp;
    
    for (i = 0; i < arrayLength - 1; i++)
    {
        for (j = i+1; j < arrayLength; j++)
        {
            if (array[i] > array[j])
            {
                temp = array[i];
                array[i] = array[j];
                array[j] = temp;
            }
        }
    }
}

