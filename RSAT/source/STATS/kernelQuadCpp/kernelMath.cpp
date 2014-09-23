#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include "kernelMath.h"
#include "sortandmerge.h"
//using namespace std;

void getXbar(double * array, int *elements, double *average){
    // calculates the average
    
    double xyVal[2];
    
    int i;
    
    average[0] = 0;
    average[1] = 0;
    
    for (i=0;i<*elements; i++) {
        
        getXY(array,&i,xyVal);
        
        average[0] = average[0] + xyVal[0];
        average[1] = average[1] + xyVal[1];
    }
    
    
    average[0] = average[0]/((double)(*elements));
    average[1] = average[1]/((double)(*elements));
    
    
    

}



void getS(double *array, int *elements, double * average, double *Sval){
    // to calcuate Standard deviation
    
    int i;
    double xyVal[2];
    
    Sval[0] = 0;
    Sval[1] = 0;
    Sval[2] = 0;
    for (i=0; i<*elements; i++) {
        getXY(array,&i,xyVal);
        
        Sval[0] = Sval[0] + (xyVal[0]-average[0])*(xyVal[0]-average[0]);
        Sval[1] = Sval[1] + (xyVal[1]-average[1])*(xyVal[1]-average[1]);
        Sval[2] = Sval[2] + (xyVal[0]-average[0])*(xyVal[1]-average[1]);
        
        
    }
    //printf("Sval %lf\n",Sval[2]);

    
    Sval[0] = Sval[0]/((double)(*elements -1));
    Sval[1] = Sval[1]/((double)(*elements -1));
    Sval[2] = Sval[2]/((double)(*elements -1));
 
    //printf("Sval %lf\n",Sval[2]);
    //printf("el %lf\n",(double)*elements);

    
    
    
    
}



void getEigenVector(double *Sval, double *EigenVector){
    double Eig[2];
    double v1,v2,norm;
    double rho,sigma1,sigma2;

    sigma1 = sqrt(Sval[0]);
    sigma2 = sqrt(Sval[1]);
    rho = Sval[2]/(sigma1*sigma2);
    
    if (fabs(rho)<1e-14){
        
        EigenVector[0] = 1;
        EigenVector[1] = 0; 
    }
    else{
        
    
    getEigenvals(Sval,Eig);

    v1 = 1;
    v2 = -(Sval[0]-Eig[0])/Sval[2];
    norm = sqrt(v1*v1 + v2*v2);
    
    EigenVector[0] = v1/norm;
    EigenVector[1] = v2/norm;
    
    
    }
    


    
    
}

void getEigenvals(double *Sval, double *Eig){
    double lam1,lam2;
    double a,b,c;
    
    
    a = 1.;
    b = -(Sval[0]+Sval[1]);
    c = Sval[0]*Sval[1]-Sval[2]*Sval[2];
         
    lam1 = (-b+sqrt(b*b-4*a*c))/(2.*a);
    lam2 = (-b-sqrt(b*b-4*a*c))/(2.*a);
    
    
    Eig[0] = lam1;
    Eig[1] = lam2;



}


void getPQ(double *array,int *elements,double *U,double *Parray,double *Qarray){
    int i;
    double P,Q;
    double xyVal[2];
    
    if (fabs(U[1])>=1e-15){
    for (i=0; i<*elements; i++) {
        getXY(array,&i,xyVal);
        P =  xyVal[0]*U[0] + xyVal[1]*U[1];
        Q = -xyVal[0]*U[1] + xyVal[1]*U[0];
        
        Parray[i]= P;
        Qarray[i] = Q;}
    }else{ 
            for (i=0; i<*elements; i++) {
                getXY(array,&i,xyVal);
                //P =  xyVal[0]*U[0] + xyVal[1]*U[1];
                //Q = -xyVal[0]*U[1] + xyVal[1]*U[0];
                
                Parray[i]= xyVal[0];
                Qarray[i] = xyVal[1];
            }
        }
        }        



void IQR(double *P, int elements, double *val){
    int mid;
    
    int low;
    int high;
    
    
    quicksortSimple(P,0,elements-1);
    
    
    
    if (elements%2==0){
        
        mid = elements/2;
        
        low = mid/2;
        
        high = elements-low;
        
        *val = .5*(P[high]+P[high+1])-.5*(P[low]+P[low+1]);

        
        
    }else{
        
        mid = (elements+1)/2;
        low = mid/2;
        high = elements - low+1;
        
        low = low-1;
        high = high -1; //-1 to work with C++ array numbering
        
        *val = P[high]-P[low];
        
    }
        
        
        
        
    
    if(*val<0){
        
        
        printf("ERROR!!!: IQR is less than zero, check kernelMath.cpp \n");
        exit(1);
        
        
    }
    
    
    
    
    
    
    
    
}


void H2parameters(double *array,int *elements,double *H2inv,double *detH2){
    
    double U[2];//eigenvector information
    double average[2];
    double Sval[3];
    double *P;
    double *Q;
    double sigma1,sigma2;
    double h1,h2;
    double IQR1,IQR2;
    
    
    P = new double[*elements];
    Q = new double[*elements];
 
    
    getXbar(array,elements, average);
    getS(array,elements,average,Sval);
    getEigenVector(Sval,U);
    getPQ(array,elements,U,P,Q);
    
    getSigma(P,elements,&sigma1);
    getSigma(Q,elements,&sigma2);
    
    IQR(P,*elements,&IQR1);
    IQR(Q,*elements,&IQR2);
    
//    printf("IQR1 %lf %lf \n",P[0],P[1]);
//    printf("IQR2 %lf %lf \n",Q[0],Q[1]);


   // printf("s1 %lf\n",sigma1);
   // printf("s2 %lf\n",sigma2);


    if (sigma1<=IQR1/1.34){
    
        h1 = 1.06*sigma1*pow((double)(*elements),-.2);
            }
    
    else{
        
        h1 = 1.06*IQR1/1.34*pow((double)(*elements),-.2);
       
    }
    
    
    if (sigma2<=IQR2/1.34){
        
        h2 = 1.06*sigma2*pow((double)(*elements),-.2);
    }
    else{
        
        h2 = 1.06*IQR2/1.34*pow((double)(*elements),-.2);
        
    }
    
    if (h1==0){
    printf("Warning h1 in KDE is zero. The points suggest that it is a non random variable (all points land in the same location \n");
    printf("This value is required for a 2x2 matrix inversion. It might return nan for desired PDF \n");
    }
    

    //printf("h1 %lf %lf %lf \n ",h1,sigma1,IQR1);
    //printf("h2 %lf %lf %lf \n",h2,sigma2,IQR2);
    *detH2 = h1*h1*h2*h2;
    
    H2inv[0] = U[0]*U[0]/(h1*h1) + U[1]*U[1]/(h2*h2);
    H2inv[1] = U[1]*U[1]/(h1*h1) + U[0]*U[0]/(h2*h2);
    H2inv[2] = U[0]*U[1]/(h1*h1)-U[0]*U[1]/(h2*h2);
    /*
    cout.precision(18);
    cout << "Mean1 \n"; 
    cout << average[0]<< endl;
    cout << average[1] << endl;
    cout << "Variance "<<endl; 
    cout << Sval[0] << endl;
    cout << Sval[1] << endl;
    cout << Sval[2] << endl;
    cout << "U"<<endl;
    cout << U[0]<<endl;
    cout << U[1] << endl;
    
    cout << "hvals" << endl;
    cout << h1 << endl; 
    cout << h2 << endl;
    
    cout << "sigma" << endl;
    cout << sigma1 << endl;
    cout << sigma2 << endl;
    cout << "mult" << endl;
    cout << pow(*elements,-.2) << endl;
    cout << IQR1 << endl;
    cout << IQR2 << endl;*/
    
}


void getSigma(double *P,int *elements,double * sigmaVal){
    
    int i;
    double xbar;
    double sigma;
    
    xbar = 0;
    sigma = 0;
    for (i=0; i<*elements; i++) {
        xbar = xbar + P[i];
    }
    
    xbar = xbar/(double)(*elements);
    
    
    for (i=0; i<*elements; i++) {
        sigma = sigma + (P[i]-xbar)*(P[i]-xbar);
    }
    
    sigma = sqrt(sigma/((double)(*elements-1)));
    *sigmaVal = sigma;
    
}



void getXY(double *array, int *index,double *xyVal){
    
    
    xyVal[0] = array[(*index)*2];
    xyVal[1] = array[(*index)*2+1];
    
}


void linspace(double begin, double end, int elements,double *result){
    double dx;
    int i;
    
    dx = (end-begin)/((double)elements-1.0);
    result[0] = begin;
    for (i=1; i<elements; i++) {
        result[i] = result[i-1]+dx; 
    }
    
    
}

