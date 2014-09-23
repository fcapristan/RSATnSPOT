#ifndef __QUADTREECME_H__
#define __QUADTREECME_H__

#include <vector>

using namespace std;

class QuadtreeCME;
class ObjectQuad;

class QuadtreeCME {
public:
    QuadtreeCME( double x, double y, double width, double height, int Nelements,double xcm = 0,double ycm =0);
    
	void					AddObject( ObjectQuad *object );
    void                    printTree();
    //void                    getVals(double, double, double);
    //vector<ObjectQuad*>      getVals(double, double, double);
   // void                     getVals(double , double , double , vector<ObjectQuad> *);
    //double                    getVals(double, double , double,double *);
    double                  getVals(double xi, double yi,double *H2inv, double q);





    
private:
	double					x;
	double					y;
	double					width;
	double					height;
    double                   xcm;
    double                   ycm;
	int						Nelements;
	vector<ObjectQuad*>			objects;
    
	QuadtreeCME *				parent;
	QuadtreeCME *				NW;
	QuadtreeCME *				NE;
	QuadtreeCME *				SW;
	QuadtreeCME *				SE;
    
    
    
	bool					Contains( QuadtreeCME *child, ObjectQuad *object );

};

#endif