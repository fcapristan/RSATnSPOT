#include <stdio.h>
#include <stdlib.h>


#include "QuadtreeCME.h"
#include "objectQuad.h"


int main(){
    
    int i;
    float x,y;
    QuadtreeCME MainTree(0,0,1000,1000,0);
   // vector<ObjectQuad> trial;
    float trial;
    float bs1;
    vector<int> trash;
    ObjectQuad localObject;//(1.1,3.4);
    ObjectQuad localObject2(1.23,3.1);

    
//    x = 3.2;
//    y = 5.1;
//    
//    localObject.x = x;
//    localObject.y =y;
     
    bs1 = 0;
    localObject.x = 1.2;
    localObject.y = 3.2;
    
    for (i=1; i<10; i++) {
        
        printf("Begin Loop\n");
        localObject.x = 1*(float)i;
        localObject.y =  -.1*(float)i;
        MainTree.AddObject(&localObject);
        
        printf("End Loop\n");
        trash.push_back(i);

        // MainTree.AddObject(&localObject2);

    }

    
   // printf("Printing Tree\n\n");
    //MainTree.printTree();
    
    printf("%d \n",(int)trash.size());     

    trial = MainTree.getVals(-10., -10., .1,&bs1);
  
    printf(" trial is %lf\n",trial);
    
    
    
        
    
    
    
    
    
    
    
    return 0;
    
}