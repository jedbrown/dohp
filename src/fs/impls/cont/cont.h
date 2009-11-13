#ifndef _CONT_H
#define _CONT_H

#include <dohpfsimpl.h>

PETSC_EXTERN_CXX_BEGIN

EXTERN dErr dFSCreate_Cont(dFS);

typedef struct {
  dTruth usecmatrix;

   /*
   * Since we don't do matrix-free computations with these elements, they don't need to be persistant (i.e. defined
   * in a structure)
   * */
#if 0 
  /* The real function space, has dofs associated with every active entity */
  dInt toff[5];     /**< Offset in unified arrays for each type */
  dInt *estart;     /**< \c estart[i] is the offset of the first degree of freedom for entity \c i */
  dInt *enodes;     /**< \c enodes[i] is the number of internal nodes associated with entity \c i */
  dInt *gstart;     /**<  Like \c estart but gives offset into global vector */
#endif

  /* The expanded function space, has dofs associated only with active regions and owned inhomogenous non-Dirichlet
  * boundary faces (i.e. those faces with non-vanishing surface integrals. */
  dInt xnents;        /**< Total number of entities in expanded space (=xtoff[dTYPE_ALL]) */
  dInt xtoff[5];      /**< Offsets for each type */
  dInt *xstart;     /**< \c xstart[xtoff[d]+i] (0<=i<xtoff[d+1]-xtoff[d]) is local index of first node of entity \c i of type \c d */

  dReal *normal;    /**< normal vector associated with non-Dirichlet boundary faces, normal[i*3+j] 0<=i<xne[2], 0<=j<3 */

  Vec x,y;          /**< vectors in expanded space */

  Mat B;
} dFS_Cont;

PETSC_EXTERN_CXX_END

#endif  /* _CONT_H */
