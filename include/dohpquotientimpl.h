#ifndef _DOHPQUOTIENTIMPL_H
#define _DOHPQUOTIENTIMPL_H

#include <dohpquotient.h>
#include <private/dmimpl.h>

// Since element quadrature and mapping contexts may have different sizes in
// different elements, we cannot simply use the PETSc data structures to handle
// local to global communication.  For now, we can just use flat arrays of
// pointers on the local process, but this will require some thought when we
// want to do dynamic load balancing.
struct _dQuotientOps {
  dErr (*update)(dQuotient);
  dErr (*setup)(dQuotient);
  dErr (*setfromoptions)(dQuotient);
  dErr (*view)(dQuotient,dViewer);
  dErr (*destroy)(dQuotient);
};


/**
* The Quotient map.  A quadrature rule and element coordinate mapping defines a quotient map on the Hilbert space
* $H^1(\Omega)$ which induces a finite topology (on the infinite dimensional space).  The same quotient map can be used
* for many different discrete function spaces $X \subset H^1$.
* */
struct p_dQuotient {
  PETSCHEADER(struct _dQuotientOps);
  dMesh                    mesh;
  dMeshTag                     qsizetag;
  dMeshESH                     loc;
  dInt                    nelems;
  dInt                   *degree;
  void                      **quad; // element quadrature context for locally owned elements
  void                      **emap; // locally owned element map context, compatible with quad
  dInt                    setupcalled;
  dQuotientSetDegreeFunc   setdegreefunc;
  void                       *setdegreectx;
  dBool                  setdegreeset;
};

extern dErr dQuotientCreate_Gauss(dQuotient);

#endif /* _DOHPQUOTIENTIMPL_H */
