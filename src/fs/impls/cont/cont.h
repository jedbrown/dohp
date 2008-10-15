#ifndef _CONT_H
#define _CONT_H

#include "private/fsimpl.h"

PETSC_EXTERN_CXX_BEGIN

extern dErr dFSCreate_Cont(dFS);

typedef struct {
  dBool usecmatrix;
  dInt D;                       /* Number of degrees of freedom per node */
  dMeshTag ruletag,degreetag;
  dInt m,n,nowned;

  /* The real function space, has dofs associated with every active entity */
  dInt ne[4];       /**< Number of entities of each type (vertices, edges, faces, regions) */
  dInt toff[4];     /**< Offset in unified arrays for each type */
  dInt *estart;     /**< \c estart[i] is the offset of the first degree of freedom for entity \c i */
  dInt *enodes;     /**< \c enodes[i] is the number of internal nodes associated with entity \c i */
  dInt *gstart;     /**<  Like \c estart but gives offset into global vector */

  /* The expanded function space, has dofs associated only with active regions and non-Dirichlet boundary faces */
  dInt xne[4];      /**< Number of expanded elements of each type, normally just inhomogenous non-Dirichlet boundary faces and regions */
  dInt xtoff[4];    /**< Offsets for each type */
  dInt *xstart;     /**< \c xstart[xtoff[d]+i] (0<=i<xne[d]) is index of first node of entity \c i of type \c d */
  dInt *xnodes;     /**< \c xnodes[xtoff[d]+i] is number of nodes associated with expanded entity \c i of type \c d */
  dRule *rule;      /**< Integration rule */
  dEFS *efs;        /**< Element function space, defined for all entities */

  dReal *normal;    /**< normal vector associated with non-Dirichlet boundary faces, normal[i*3+j] 0<=i<xne[2], 0<=j<3 */

  Sliced sliced;    /**< Manages the relationship between global vector and local vector */
  Vec x,y;          /**< vectors in expanded space */

  Mat F;              /**< facet constraint matrix */
  Mat Fp;             /**< facet preconditioning constraint matrix */
  Mat C;              /**< full-order constraint matrix (element dofs to local numbering) */
  Mat Cp;             /**< preconditioning constraint matrix (element dofs to local numbering, as sparse as possible) */
} dFS_Cont;

PETSC_EXTERN_CXX_END

#endif  /* _CONT_H */
