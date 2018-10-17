#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define W 4

typedef struct komvos* deiktis;
typedef struct komvos
{
	deiktis next;
    int indexNo;
    long FiValue;
}komvos;

void katastrofi(deiktis *arxi);
void eisagogi_arxi(deiktis *arxi, int thesi);

