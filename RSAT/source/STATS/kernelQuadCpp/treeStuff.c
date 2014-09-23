#include <stdio.h>
#include <stdlib.h>


 struct node{
    double x;
    double y;
    int nElements;
    struct node *child[4]; //
    // array order NE -> 0, NW -> 1, SW -> 2, SE ->3 
 } ;

void InitializeTree(tree *);





int main(){
    int i;
    tree quadtree; 
    //quadtree = new tree;
    tree *temp;
    
   // quadtree = new tree;
    InitializeTree(&quadtree);
    
    printf("%d\n",quadtree.nElements);


    
    
    
    return 0;
    
}


struct node* NewNode(double x, double y){
    
    struct node* node = new(struct node);
    node->x = x;
    node->y = y;
    node->child[0]=NULL;
    node->child[1]=NULL;
    node->child[2]=NULL;
    node->child[3]=NULL;

    return node
    
}







