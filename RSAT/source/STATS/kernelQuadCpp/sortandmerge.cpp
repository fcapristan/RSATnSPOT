#include <math.h>
#include "sortandmerge.h"

void quicksort(double *array,int leftIndex, int rightIndex){
    
    int pivotIndex;
    int pivotNewIndex;
    
    
    
    if (leftIndex<rightIndex){
        
        
        
        pivotIndex = floor(.5*((float)leftIndex+(float)rightIndex));
        
        
        partition(array,leftIndex,rightIndex,pivotIndex,&pivotNewIndex);
        quicksort(array,leftIndex,pivotNewIndex-1);
        quicksort(array,pivotNewIndex+1,rightIndex);
        
        
    }
    
    
}

void quicksortSimple(double *array,int leftIndex, int rightIndex){
    
    int pivotNewIndex;
    int pivotIndex;
    
    
    
    if (leftIndex<rightIndex){
        
        
        
        pivotIndex = floor(.5*((float)leftIndex+(float)rightIndex));
        
        
        partitionSimple(array,leftIndex,rightIndex,pivotIndex,&pivotNewIndex);
        quicksortSimple(array,leftIndex,pivotNewIndex-1);
        quicksortSimple(array,pivotNewIndex+1,rightIndex);
        
        
    }
    
    
    
    
    
    
    
    
}



void partitionSimple(double *array,int leftIndex, int rightIndex, int pivotIndex,int *storeIndex){
    int i;
    double pivotValue,localValue;
    double temp;
    pivotValue = array[pivotIndex];
    //valsFromIndex(array,&pivotIndex,&pivotValue);
    //swapArray(array,&pivotIndex,&rightIndex);
    temp = array[pivotIndex];
    array[pivotIndex] = array[rightIndex];
    array[rightIndex] = temp;
    
    *storeIndex = leftIndex;
    
    for (i=leftIndex;i<rightIndex;i++){
        
        //valsFromIndex(array,&i,&localValue);
        localValue = array[i];
        
        if (localValue<pivotValue){
            
            
            
            //swapArray(array,&i,storeIndex);
            
            temp = array[i];
            array[i] = array[*storeIndex];
            array[*storeIndex] = temp;
            *storeIndex = *storeIndex + 1;
            
            
        }
        
        
        
        
        
    }
    
    
    
    //swapArray(array,storeIndex,&rightIndex);
    temp = array[*storeIndex];
    array[*storeIndex] = array[rightIndex];
    array[rightIndex] = temp;
    
    
}
















void partition(double *array,int leftIndex, int rightIndex, int pivotIndex,int *storeIndex){
    int i;
    double pivotValue,localValue;
    
    
    valsFromIndex(array,&pivotIndex,&pivotValue);
    swapArray(array,&pivotIndex,&rightIndex);
    *storeIndex = leftIndex;
    
    for (i=leftIndex;i<rightIndex;i++){
        
        valsFromIndex(array,&i,&localValue);
        
        if (localValue<pivotValue){
            
            
            
            swapArray(array,&i,storeIndex);
            *storeIndex = *storeIndex + 1;
            
            
        }
        
        
        
        
        
    }
    
    
    
    swapArray(array,storeIndex,&rightIndex);
    
 
    
}


void valsFromIndex(double *array,int *index,double *val){
    
    // gets distance value for the desired index (xy location) 
    
    *val = array[(*index)*3+2];
    
    
    
}


void assignValsFromIndex(double *a1, int *index1,double *result, int *index2){
    
    result[(*index2)*3] = a1[(*index1)*3];
    result[(*index2)*3+1] = a1[(*index1)*3+1];
    result[(*index2)*3+2] = a1[(*index1)*3+2];
    
    
}



void swapArray(double *array,int *index1,int *index2){
    
    double temp[3];
    
    
    
    temp[0] = array[(*index1)*3];
    temp[1] = array[(*index1)*3 + 1];
    temp[2] = array[(*index1)*3 + 2];
    
    
    array[(*index1)*3] = array[(*index2)*3];
    array[(*index1)*3 + 1] = array[(*index2)*3 + 1];
    array[(*index1)*3 + 2] = array[(*index2)*3 + 2];
    
    
    
    array[(*index2)*3] = temp[0];
    array[(*index2)*3 + 1] = temp[1];
    array[(*index2)*3 + 2] = temp[2];
    
    
    
    
    
}



void mergedata(double *array1,int len1, double *array2,int len2,double *result){
    int i,j,k;
    double val1,val2;
    
    i =0;j=0;k=0;
    
    
    
    
    while(i<len1 && j<len2){
        
        
        
        valsFromIndex(array1,&i,&val1);
        valsFromIndex(array2,&j,&val2);
        
        if(val1<val2){
            
            
            
            assignValsFromIndex(array1, &i,result, &k);
            
            
            //result[k] = array1[i];
            i++;k++;
            
        }
        else{
            
            
            assignValsFromIndex(array2, &j,result, &k);
            
            //    result[k] = array2[j]
            j++;k++;
            
        }
    }
    
    
    if(i==len1 && j<len2){
        
        while(j<len2){
            
            
            assignValsFromIndex(array2, &j,result, &k);
            
            //                result[k] = array2[j];
            
            j++;k++;
            
            
        }
        
        
    }
    else if(j==len2 && i<len1){
        while(i<len1){
            //       printf("k is %d \n",k);
            //      printf("sum %d \n",len1+len2);
            
            assignValsFromIndex(array1, &i,result, &k);
            
            
            
            //                result[k] = array1[i];
            i++;k++;
            
            
        }
        
        
        
    }
    
    
    
    
}
