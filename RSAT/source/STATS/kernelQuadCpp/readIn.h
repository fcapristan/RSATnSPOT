#ifndef READIN_H_
#define READIN_H_

#include <stdio.h>
#include <stdlib.h>

void reader(FILE * ifp,int *startloc,int *endloc, double *xydist);
void getStartlocEndloc(int rank,int size, int lines,int *startloc,int *endloc);
void getElements(int rank,int size,int lines,int *elements);


#endif