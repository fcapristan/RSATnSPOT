#ifndef SERIALKERNEL_H
#define SERIALKERNEL_H

int KernelQuadTree(int N, double *xy,int rows,int cols,double *Xmesh,double *Ymesh,double *fMesh,double q,int QuadTreeCalc);
int weightedKDE(int N, double *xy,int rows,int cols,double *Xmesh,double *Ymesh,double *fMesh,double *weight);



#endif
