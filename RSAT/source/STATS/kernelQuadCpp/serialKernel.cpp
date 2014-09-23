#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kernelMath.h"
#include "sortandmerge.h"
#include "quadtreeCME.h"
#include "objectQuad.h"
#include "serialKernel.h"
using namespace std;


int KernelQuadTree(int N, double *xy,int rows,int cols,double *Xmesh,double *Ymesh,double *fMesh,double q,int QuadTreeCalc){
 // N corresponds to the number of samples in xy 
    // xy is longitude and latitude (Nx2 Matrix)
    // Xmesh and Ymesh -> rowsxcols matrix
    // q -> parameter required for quadtree implementation q = D/r
    // QuadTreeCalc -> 0 for regular KDE, 1 for Quadtree implementation
    const double PI  =3.141592653589793238462643;
    // values needed for H2 determination
    double H2inv[3];
    double detH2,xHx;
    // index for for loops
    int i,j,k,step;
    double val,xi,yi,xcm,ycm; // temporary value to hold a sum
    double denom;// just to avoid mutiple calculations
    QuadtreeCME MainTree(0,0,600,600,0); // starting tree centered at 0 with 600 height and width...height and width MUST BE EQUAL
    ObjectQuad localObject;
    
    
    H2parameters(xy,&N,H2inv,&detH2);
    //printf("detH2 %lf\n",detH2);
    //printf("H2 %lf\n",H2inv[0]);
    //printf("H2 %lf\n",H2inv[1]);
    //printf("H2 %lf\n",H2inv[2]);

    denom = 2.0*PI*sqrt(detH2)*(double)N;


    //printf("denom %lf\n",denom);
    if(QuadTreeCalc==0){
        for (i=0;i<rows;i++){
            for (j=0;j<cols;j++){
                xi = Xmesh[cols*i+j];
                yi = Ymesh[cols*i+j];
                val = 0; // resetting val to zero for each grid
                for (k=0;k<N;k++){
                    xcm = xy[2*k];
                    ycm = xy[2*k+1];
                    xHx = (xi-xcm)*(xi-xcm)*H2inv[0] + 2*(xi-xcm)*(yi-ycm)*H2inv[2] + (yi-ycm)*(yi-ycm)*H2inv[1];
                    
                    val = exp(-.5*xHx)+val;
                    
                }
                fMesh[cols*i+j] = val/denom;
                //printf("fvals %lf\n",val);

                
            }
                                    
        }
    }
    else if(QuadTreeCalc==1){
        
        for (i=0;i<N;i++){
            localObject.x = xy[2*i];
            localObject.y = xy[2*i+1];
            MainTree.AddObject(&localObject);         
        }
        
        for (i=0;i<rows;i++){
            for (j=0;j<cols;j++){
                xi = Xmesh[cols*i+j];
                yi = Ymesh[cols*i+j];
                val = MainTree.getVals(xi,yi,H2inv,q);
                fMesh[cols*i+j]=val/denom;
               // printf("fvals %lf\n",val);

  
            }
        
        }
        
    }
}




int weightedKDE(int N, double *xy,int rows,int cols,double *Xmesh,double *Ymesh,double *fMesh,double *weight){
    // N corresponds to the number of samples in xy 
    // xy is longitude and latitude (Nx2 Matrix)
    // Xmesh and Ymesh -> rowsxcols matrix
    // q -> parameter required for quadtree implementation q = D/r
    // QuadTreeCalc -> 0 for regular KDE, 1 for Quadtree implementation
    const double PI  =3.141592653589793238462643;
    // values needed for H2 determination
    double H2inv[3];
    double detH2,xHx;
    // index for for loops
    int i,j,k,step;
    double val,xi,yi,xcm,ycm; // temporary value to hold a sum
    double denom;// just to avoid mutiple calculations
    QuadtreeCME MainTree(0,0,600,600,0); // starting tree centered at 0 with 600 height and width...height and width MUST BE EQUAL
    ObjectQuad localObject;
    
    
    H2parameters(xy,&N,H2inv,&detH2);
    //printf("detH2 %lf\n",detH2);
    //printf("H2 %lf\n",H2inv[0]);
    //printf("H2 %lf\n",H2inv[1]);
    //printf("H2 %lf\n",H2inv[2]);
    
    denom = 2.0*PI*sqrt(detH2)*(double)N;
    
    
        for (i=0;i<rows;i++){
            for (j=0;j<cols;j++){
                xi = Xmesh[cols*i+j];
                yi = Ymesh[cols*i+j];
                val = 0; // resetting val to zero for each grid
                for (k=0;k<N;k++){
                    xcm = xy[2*k];
                    ycm = xy[2*k+1];
                    xHx = (xi-xcm)*(xi-xcm)*H2inv[0] + 2*(xi-xcm)*(yi-ycm)*H2inv[2] + (yi-ycm)*(yi-ycm)*H2inv[1];
                    
                    val = weight[k]*exp(-.5*xHx)+val;
                    
                }
                fMesh[cols*i+j] = val/denom;
                //printf("fvals %lf\n",val);
                
                
            }
            
        }
    

}
