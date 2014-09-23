# distutils: language = c++
# distutils: sources = serialKernel.cpp kernelMath.cpp sortandMerge.cpp quadtreeCME.cpp objectQuad.cpp

import numpy as np
cimport numpy as np

cdef extern from "serialKernel.h" :
     int KernelQuadTree(int N, double *xy,int rows,int cols,double *Xmesh,double *Ymesh,double *fMesh,double q,int QuadTreeCalc)
     int weightedKDE(int N, double *xy,int rows,int cols,double *Xmesh,double *Ymesh,double *fMesh,double *weight)

def serialK(N,xy,rows,cols,Xmesh,Ymesh,q,QuadtreeCalc):
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] xy_c
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] Xmesh_c
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] Ymesh_c
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] fMesh_c
    xy_c = np.ascontiguousarray(xy, dtype=np.float64)
    Xmesh_c = np.ascontiguousarray(Xmesh, dtype=np.float)
    Ymesh_c = np.ascontiguousarray(Ymesh, dtype=np.float)
    fMesh = np.zeros(([rows,cols]))
    fMesh_c = np.ascontiguousarray(fMesh, dtype=np.float)
    KernelQuadTree( N,&xy_c[0,0],rows,cols,&Xmesh_c[0,0],&Ymesh_c[0,0],&fMesh_c[0,0],q,QuadtreeCalc)
    return fMesh_c

def serialWeightKDE(N,xy,rows,cols,Xmesh,Ymesh,weights):
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] xy_c
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] Xmesh_c
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] Ymesh_c
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] fMesh_c
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] weights_c
    weights_c = np.ascontiguousarray(weights,dtype=np.float)
    xy_c = np.ascontiguousarray(xy, dtype=np.float64)
    Xmesh_c = np.ascontiguousarray(Xmesh, dtype=np.float)
    Ymesh_c = np.ascontiguousarray(Ymesh, dtype=np.float)
    fMesh = np.zeros(([rows,cols]))
    fMesh_c = np.ascontiguousarray(fMesh, dtype=np.float)
    weightedKDE( N,&xy_c[0,0],rows,cols,&Xmesh_c[0,0],&Ymesh_c[0,0],&fMesh_c[0,0],&weights_c[0])
    return fMesh_c

