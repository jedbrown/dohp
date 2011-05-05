#!/usr/bin/env python2
from __future__ import division

import sympy
from sympy import Matrix, ccode
from dohpexact import *

def SecondInvariant(Du):
    return 0.5*Du.dot(Du)

class StokesExact(Exact):
    def __init__(self, name=None, model='B eps pe', param='a b c'):
        Exact.__init__(self, name=name, model=model, param=param, fieldspec=[('u',3), ('p',1)])
    def eta(self, gamma):       # Power law
        B, eps, pe = self.model_get('B eps pe')
        return B * (0.5*eps**2 + gamma)**((pe-2)/2)
    def weak_homogeneous(self, x, U, dU, V, dV):
        u, p = U[:3], U[3]
        Du = sym(dU[:3,:])
        v, q = V[:3], V[3]
        Dv = sym(dV[:3,:])
        gamma = SecondInvariant(Du)
        eta = self.eta(gamma)
        return eta*Dv.dot(Du) - q*Du.trace() - p*Dv.trace()
    def create_prototype(self, definition=False):
        return 'dErr StokesCaseCreate_%(name)s(StokesCase case)%(term)s' % dict(name=self.name, term=('' if definition else ';'))
    def solution_prototype(self):
        return 'static dErr StokesCaseSolution_%(name)s(StokesCase scase,const dReal x[3],dScalar u[3],dScalar du[9],dScalar p[1],dScalar dp[3])' % dict(name=self.name)
    def forcing_prototype(self):
        return 'static dErr StokesCaseForcing_%(name)s(StokesCase scase,const dReal x[3],dScalar fu[3],dScalar fp[1])' % dict(name=self.name)
    def solution_code(self):
        from sympy.abc import a,b,c
        x = Matrix(symbol3('x'))
        def body():
            for (i,expr) in enumerate(self.solution_scaled(x,a,b,c)):
                yield ccode(expr, assign_to=self.fieldname(i))
            for (i,expr) in enumerate(self.solution_gradient(x,a,b,c)):
                yield ccode(expr, assign_to=self.dfieldname(i))
        return '''
%(prototype)s
{
  const StokesCase_Exact *ctx = scase->data;
  const dReal a = ctx->a,b = ctx->b,c = ctx->c,scale = ctx->scale;
  %(body)s
  return 0;
}''' % dict(prototype=self.solution_prototype(), body='\n  '.join(body()))
    def forcing_code(self):
        from sympy.abc import a,b,c
        x = Matrix(symbol3('x'))
        def body():
            U = self.solution_scaled(x,a,b,c)
            dU = self.solution_gradient(x,a,b,c)
            V,dV = self.residual(x,U,dU)
            for i,row in enumerate(rows(dV)):
                yield ccode(V[i] - divergence(row,x), assign_to=self.ffieldname(i))
        return '''
%(prototype)s
{
  const StokesCase_Exact *ctx = scase->data;
  const dReal a = ctx->a,b = ctx->b,c = ctx->c,scale = ctx->scale;
  const dReal B = scase->rheo.B,eps = scase->rheo.eps,pe = scase->rheo.p;
  %(body)s
  return 0;
}
''' % dict(prototype=self.forcing_prototype(), body='\n  '.join(body()))
    def create_code(self):
        return '''
static dErr StokesCaseCreate_%(name)s(StokesCase scase)
{
  dErr err;
  StokesCase_Exact *ctx;

  dFunctionBegin;
  err = dCalloc(sizeof(*ctx),&ctx);dCHK(err);
  ctx->a     = 1.0;
  ctx->b     = 1.0;
  ctx->c     = 1.0;
  ctx->scale = 1.0;
  scase->data = ctx;
  scase->solution       = StokesCaseSolution_%(name)s;
  scase->forcing        = StokesCaseForcing_%(name)s;
  scase->setfromoptions = StokesCaseSetFromOptions_Exact;
  scase->destroy        = StokesCaseDestroy_Exact;
  dFunctionReturn(0);
}
''' % dict(name=self.name)

class Exact0(StokesExact):
    'Classical 2D converging flow'
    def solution(self, x,y,z, a,b,c):
        from sympy import sin, cos, pi
        return Matrix([ a * sin(pi*x/2) * cos(pi*y/2),
                       -b * cos(pi*x/2) * sin(pi*y/2),
                        1 * (c-1) * z,
                        0.25*(cos(pi*x) + cos(pi*y)) + 10*(x+y)])
class Exact1(StokesExact):
    'From Larin & Reusken, 2009, with pressure shifted by a constant'
    def solution(self, x,y,z, a,b,c):
        from sympy import sin,cos,pi
        xx, yy, zz = (pi*s/2 for s in [x,y,z])
        return Matrix([+a/3 * sin(xx) * sin(yy) * sin(zz),
                       -b/3 * cos(xx) * cos(yy) * sin(zz),
                       -c*2/3 * cos(xx) * sin(yy) * cos(zz),
                        1 + cos(xx) * sin(yy) * sin(zz)])
class Exact2(StokesExact):
    def solution(self, x,y,z, a,b,c):
        return Matrix([a*z**3,
                       b*x**3,
                       c*y**3,
                       (1-x**2)*(1-y**2)*(1-z**2)])
class Exact3(StokesExact):
    def solution(self, x,y,z, a,b,c):
        from sympy import sin,cos,pi
        xx, yy, zz = (pi*s/2 for s in [x,y,z])
        return Matrix([a*sin(2*xx),
                       b*cos(yy),
                       c*cos(zz),
                       cos(xx)*cos(yy)*cos(zz)])

def implementation():
  return '''
#include <stokesimpl.h>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

typedef struct {
  dReal a,b,c,scale;
} StokesCase_Exact;

static dErr StokesCaseSetFromOptions_Exact(StokesCase scase) {
  StokesCase_Exact *exc = scase->data;
  dErr err;

  dFunctionBegin;
  err = PetscOptionsHead("StokesCase_Exact options");dCHK(err); {
    err = PetscOptionsReal("-exact_a","First scale parameter","",exc->a,&exc->a,NULL);dCHK(err);
    err = PetscOptionsReal("-exact_b","Second scale parameter","",exc->b,&exc->b,NULL);dCHK(err);
    err = PetscOptionsReal("-exact_c","Third scale parameter","",exc->c,&exc->c,NULL);dCHK(err);
    err = PetscOptionsReal("-exact_scale","Overall scale parameter","",exc->scale,&exc->scale,NULL);dCHK(err);
  } err = PetscOptionsTail();dCHK(err);
  dFunctionReturn(0);
}
static dErr StokesCaseDestroy_Exact(StokesCase scase) {
  dErr err;

  dFunctionBegin;
  err = dFree(scase->data);dCHK(err);
  dFunctionReturn(0);
}
'''

def registerall(list):
    def register(sol):
        return 'err = StokesCaseRegister("%(name)s",StokesCaseCreate_%(name)s);dCHK(err);' % dict(name=sol.name)
    return '''
dErr StokesCaseRegisterAll_Exact(void)
{
  dErr err;
  dFunctionBegin;
  %s
  dFunctionReturn(0);
}
''' % '\n  '.join((register(sol) for sol in list))

def exact_body(name):
    return '''

'''

if __name__ == "__main__":
    import pdb
    solutions = [Exact0(), Exact1(), Exact2(), Exact3()]
    with open('stokesexact.c', 'w') as fimpl:
        fimpl.write(implementation())
        for sol in solutions:
            fimpl.write(sol.solution_code())
            fimpl.write(sol.forcing_code())
            fimpl.write(sol.create_code())
        fimpl.write(registerall(solutions))
    print('Wrote stokesexact.c')
