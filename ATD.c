#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ATD.h"

void katastrofi(deiktis *arxi)
{
    deiktis temp;
    
    while (*arxi != NULL)
    {
        temp = *arxi;
        *arxi = temp->next;
        free(temp);
    }
    *arxi = NULL;
}

void eisagogi_arxi(deiktis *arxi, int thesi)
{
    deiktis neos = malloc(sizeof(komvos));
    if (neos == NULL)
    {
        printf("...\n");
        exit(-1);
    }
    
   	neos->indexNo = thesi;
    neos->next = *arxi;
    *arxi = neos;
}
