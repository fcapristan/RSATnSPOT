#include <math.h>
#include "readIn.h"

using namespace std;




void reader(FILE * ifp,int *startloc,int *endloc, double *xydist){
    
    // this function reads in the xy location file. It also calculates and stores the distance
    
    // with respect to (zero,zero). It returns am array of the form 'x,y,distance'
    
    double val1,val2;//this will temporarily store the values x y
    int counter,countlines;
    int i;
    
    
    
    counter = 0;
    
    for (i=0;i<=*endloc;i++){
        
        fscanf(ifp,"%lf %lf",&val1,&val2);//read in values from x y location file
        
        
        if (i>=*startloc){
            
            xydist[counter] = val1;            
            xydist[counter+1] = val2;            
            xydist[counter+2] = sqrt(val1*val1+val2*val2);

            
            counter = counter+3;
            
            
            
        }
        
        
        
    }
    
    
    
}







void getStartlocEndloc(int rank,int size, int lines,int *startloc,int *endloc){
    
    int nelements;
    int nleft;
    
    
    nelements = int(floor(float(lines)/float(size)));
    nleft = (int)lines%size;
    
    
    
    if (rank<nleft){
        
        nelements = nelements+1;
        *startloc = rank*(nelements);
        
        
    }
    else{
        
        *startloc = rank*(nelements)+nleft;
        
        
    }
    
    // printf("left %d %d \n",nelements,size);
    
    *endloc = *startloc + nelements -1 ;
    
}






void getElements(int rank,int size,int lines,int *elements){
    
    
    int nelements;
    int nleft;
    
    
    nelements = int(floor(float(lines)/float(size)));
    nleft = (int)lines%size;
    
    
    
    if (rank<nleft){
        
        nelements = nelements+1;
        
        
    }
    
    *elements = nelements;  
    
    
    
    
    
}