#include <math.h>
#include <stdio.h>
#include "quadtreeCME.h"
#include "objectQuad.h"
using namespace std;

QuadtreeCME::QuadtreeCME( double _x, double _y, double _width, double _height, int _Nelements, double _xcm,double _ycm) :
x		( _x ),
y		( _y ),
width	( _width ),
height	( _height ),
Nelements	( _Nelements ),
xcm         (_xcm),
ycm         (_ycm)
{
//	if ( Nelements <2 ) {
//		return; 
//	}
//    
//	NE = new Quadtree( x + width/4.0f, y+width/4.0f, width / 2.0f, height / 2.0f, 0 );
//    NW = new Quadtree( x - width/4.0f, y+width/4.0f, width / 2.0f, height / 2.0f, 0 );
//    
//    SE = new Quadtree( x + width/4.0f, y-width/4.0f, width / 2.0f, height / 2.0f, 0 );
//    SW = new Quadtree( x - width/4.0f, y-width/4.0f, width / 2.0f, height / 2.0f, 0 );

}


void QuadtreeCME::AddObject(ObjectQuad *object ) {
    
    
    

    
	if ( Nelements==0 ) {
        
        
		objects.push_back( object );
        Nelements = 1;
        
        xcm = object->x;
        ycm = object->y;
        
    //    printf("%f %f %d \n",xcm,ycm,Nelements);

		return;
	}
    if(Nelements==1){
        ObjectQuad temp;
        
        temp.x = xcm;
        temp.y = ycm;
        
        NE = new QuadtreeCME( x + width/4.0, y+height/4.0, width / 2.0, height / 2.0, 0 );
        NW = new QuadtreeCME( x - width/4.0, y+height/4.0, width / 2.0, height / 2.0, 0 );
        
        SE = new QuadtreeCME( x + width/4.0, y-height/4.0, width / 2.0, height / 2.0, 0 );
        SW = new QuadtreeCME( x - width/4.0, y-height/4.0, width / 2.0, height / 2.0, 0 );
        
        
        
        if((object->x==xcm)&&(object->y==ycm)){
            
            printf("ERROR : 2 values cannot be the same\n");
            
            
        }
     
        
        
        
        
        
        
        if ( Contains( NW, &temp ) ) {
            NW->AddObject( &temp ); 
            NW->Nelements = 1;
            // return;
        } else if ( Contains( NE, &temp ) ) {
            NE->AddObject( &temp ); 
            NE->Nelements = 1;

            //   return;
        } else if ( Contains( SW, &temp ) ) {
            SW->AddObject( &temp ); 
            SW->Nelements = 1;

            //   return;
        } else if ( Contains( SE, &temp ) ) {
            SE->AddObject( &temp );
            SE->Nelements = 1;

            //   return;
        }else{
            
            printf("ERROR:Out of Bounds, check domain\n");
            
            
        }
        
        
        

        
        
        
        
        
        
        
    }
        

   // printf("whee\n");
   // printf("%f \n",object->x);
    
    Nelements = Nelements+1;

    
    xcm = (xcm*(double)(Nelements-1) + object->x)/((double)Nelements);
    ycm = (ycm*(double)(Nelements-1) + object->y)/((double)Nelements);
    
  //  printf("%f %f %d \n",xcm,ycm,Nelements);


    
    
    
    
    
	if ( Contains( NW, object ) ) {
		NW->AddObject( object ); 
        return;
	} else if ( Contains( NE, object ) ) {
		NE->AddObject( object ); 
        return;
	} else if ( Contains( SW, object ) ) {
		SW->AddObject( object ); 
       return;
	} else if ( Contains( SE, object ) ) {
		SE->AddObject( object );
        return;
	}else{
        
        printf("ERROR:Out of Bounds, check domain\n");
        
        
    }
    
    
    
    
    
//	if ( Contains( this, object ) ) {
//		objects.push_back( object );
//	}
}



void QuadtreeCME::printTree(){
    
    if (Nelements==0){
        //ignore
       // printf("Empty leaf\n");
        
        return;
        
    }
    if (Nelements==1){
        printf("%f %f %d \n",xcm,ycm,Nelements);
        return;
    }
    
    printf("%f %f %d \n",xcm,ycm,Nelements);
    NW->printTree();
    NE->printTree();
    SW->printTree();
    SE->printTree();

    
    
}






double  QuadtreeCME::getVals(double xi, double yi,double *H2inv, double q){
    
    double r;
    double val;
    double xHx;
    //vector<ObjectQuad*> returnObjects;
    ObjectQuad childReturnQuad;
    
    
    childReturnQuad.x = xcm;
    childReturnQuad.y = ycm;
    childReturnQuad.Nelements = Nelements;
    
    if (Nelements==0){
        //ignore
        // printf("Empty leaf\n");
        
        return 0;
        
    }
    
    
    
    if (Nelements==1){
      //printf("cond 1 left %f %f %d \n",xcm,ycm,Nelements);
        
        //returnObjects.push_back(childReturnQuad);
       // printf("width  %lf \n",width);
        
        xHx = (xi-xcm)*(xi-xcm)*H2inv[0] + 2*(xi-xcm)*(yi-ycm)*H2inv[2] + (yi-ycm)*(yi-ycm)*H2inv[1];
        val = Nelements*(exp(-.5*xHx));
        return val;// returnObjects;
        
    }
    
    r = sqrt(pow(xcm-xi,2)+pow(ycm-yi,2));

  //  printf("r  %lf \n",r);

    
    if(width/r<q){
        
      //  printf("cond %f %f %d \n",xcm,ycm,Nelements);
        
      //  returnObjects.push_back(childReturnQuad);
        //val = xcm+ ycm+(*param);
        
        
        xHx = (xi-xcm)*(xi-xcm)*H2inv[0] + 2*(xi-xcm)*(yi-ycm)*H2inv[2] + (yi-ycm)*(yi-ycm)*H2inv[1];
        val = Nelements*(exp(-.5*xHx));
        

        return val;// returnObjects;
                             
                             
        }
    
    
    
    
                           
    val = NW->getVals(xi,yi,H2inv,q)+
    NE->getVals(xi,yi,H2inv,q)+
    SW->getVals(xi,yi,H2inv,q)+
    SE->getVals(xi,yi,H2inv,q);    
    
    return val;
                             
}
                             
                             




bool QuadtreeCME::Contains( QuadtreeCME *child, ObjectQuad *object ) {
	
    double lowx,lowy,highx,highy;
    
    
    lowx = child->x - child->width/2.0;
    lowy = child->y - child->height/2.0;
    highx = child->x + child->width/2.0;
    highy = child->y + child->height/2.0;
    
   // printf("lims %f %f %f %f %f %f\n",lowx,highx,lowy,highy,object->x,object->y);
    
    return	 ( object->x <= highx &&
               object->y <= highy &&
               object->x >= lowx  &&
               object->y >= lowy);
    
    
        
    
    
    
    
}
