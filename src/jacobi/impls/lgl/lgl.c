/**
* @file   lgl.c
* @author Jed Brown <jed@59A2.org>
* @date   Fri Aug 22 20:39:12 2008
* 
* @brief  A nodal Legendre-Gauss-Lobatto basis and Gauss quadrature.
* 
* 
*/

#include "petsc.h"
#include "dohpjacobi.h"

#undef __FUNCT__
#define __FUNCT__ "DohpJacobiCreate_LGL"
PetscErrorCode DohpJacobiCreate_LGL(DohpJacobi jac)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscPrintf(((PetscObject)jac)->comm,"DohpJacobiCreate_LGL (need to do something here)\n");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

