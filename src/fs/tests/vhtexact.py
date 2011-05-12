#!/usr/bin/env python2
from __future__ import division

import sympy
from sympy import Matrix, ccode
from dohpexact import *

def SecondInvariant(Du):
    return 0.5*Du.dot(Du)
def transition(a, b, width, x): # Smooth function transitioning from a to b near 0 with characteristic width
    from sympy import tanh
    return a + (b-a)*((1+tanh(x/width))/2)

class VHTExact(Exact):
    def __init__(self, name=None, model='B0 R Q eps pe kappa0 kappa1 T0', param='a b c'):
        Exact.__init__(self, name=name, model=model, param=param, fieldspec=[('u',3), ('p',1), ('e',1)])
    def eta(self, gamma, e):       # Power law with Arrhenius relation
        from sympy import exp
        B0, R, Q, eps, pe = self.model_get('B0 R Q eps pe')
        n = 1/(pe - 1)
        T = self.temperature(e)
        B = B0 * exp(Q / (n*R*T))
        return B * (0.5*eps**2 + gamma)**((pe-2)/2)
    def temperature(self, e):
        T0, = self.model_get('T0')
        return T0 + e
    def kappa(self, e):
        kappa0,kappa1 = self.model_get('kappa0 kappa1')
        return transition(kappa0, kappa1, 1, e)
    def weak_homogeneous(self, x, U, dU, V, dV):
        u, p, e = U[:3,:], U[3], U[4]
        Du = sym(dU[:3,:])
        de = dU[4,:].T
        u_, p_, e_ = V[:3,:], V[3], V[4]
        Du_ = sym(dV[:3,:])
        de_ = dV[4,:].T
        gamma = SecondInvariant(Du)
        eta = self.eta(gamma, e)
        kappa = self.kappa(e)
        heatflux = -kappa * de + e * u
        Sigma = Du.dot(eta*Du)
        stokes = Du_.dot(eta * Du) - p_*Du.trace() - p*Du_.trace()
        enthalpy = e_ * (-Sigma) - de_.dot(heatflux)
        return stokes + enthalpy
    def create_prototype(self, definition=False):
        return 'dErr VHTCaseCreate_%(name)s(VHTCase case)%(term)s' % dict(name=self.name, term=('' if definition else ';'))
    def solution_prototype(self):
        return 'static dErr VHTCaseSolution_%(name)s(VHTCase scase,const dReal x[3],dScalar u[3],dScalar du[9],dScalar p[1],dScalar dp[3],dScalar e[1],dScalar de[1])' % dict(name=self.name)
    def forcing_prototype(self):
        return 'static dErr VHTCaseForcing_%(name)s(VHTCase scase,const dReal x[3],dScalar fu[3],dScalar fp[1],dScalar fe[1])' % dict(name=self.name)
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
  const VHTCase_Exact *ctx = scase->data;
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
  const VHTCase_Exact *ctx = scase->data;
  const dReal %(param)s;
  const dUNUSED dReal %(rheo)s;
  %(body)s
  return 0;
}
''' % dict(prototype=self.forcing_prototype(),
           param=','.join('%s = ctx->%s' % (v,v) for v in self._param.keys()),
           rheo=','.join('%s = scase->rheo.%s' % (v,v) for v in self._model.keys()),
           body='\n  '.join(body()))
    def create_code(self):
        return '''
static dErr VHTCaseCreate_%(name)s(VHTCase scase)
{
  dErr err;
  VHTCase_Exact *ctx;

  dFunctionBegin;
  err = dCalloc(sizeof(*ctx),&ctx);dCHK(err);
  %(param)s;
  scase->data = ctx;
  scase->solution       = VHTCaseSolution_%(name)s;
  scase->forcing        = VHTCaseForcing_%(name)s;
  scase->setfromoptions = VHTCaseSetFromOptions_Exact;
  scase->destroy        = VHTCaseDestroy_Exact;
  dFunctionReturn(0);
}
''' % dict(name=self.name,
           param=';\n  '.join('ctx->%s = 1.0' % v for v in self._param.keys()))

class Exact0(VHTExact):
    'Classical 2D converging flow'
    def solution(self, x,y,z, a,b,c):
        from sympy import sin, cos, pi
        return Matrix([ a * sin(pi*x/2) * cos(pi*y/2),
                       -b * cos(pi*x/2) * sin(pi*y/2),
                        1 * (c-1) * z,
                        0.25*(cos(pi*x) + cos(pi*y)) + 10*(x+y),
                        sin(pi*y/2) * cos(pi*z/2)])
class Exact1(VHTExact):
    'From Larin & Reusken, 2009, with pressure shifted by a constant'
    def solution(self, x,y,z, a,b,c):
        from sympy import sin,cos,pi
        xx, yy, zz = (pi*s/2 for s in [x,y,z])
        return Matrix([+a/3 * sin(xx) * sin(yy) * sin(zz),
                       -b/3 * cos(xx) * cos(yy) * sin(zz),
                       -c*2/3 * cos(xx) * sin(yy) * cos(zz),
                        1 + cos(xx) * sin(yy) * sin(zz),
                        sin(pi*y/2) * cos(pi*z/2)])
class Exact2(VHTExact):
    def solution(self, x,y,z, a,b,c):
        return Matrix([a*z**3,
                       b*x**3,
                       c*y**3,
                       (1-x**2)*(1-y**2)*(1-z**2)])
class Exact3(VHTExact):
    def solution(self, x,y,z, a,b,c):
        from sympy import sin,cos,pi
        xx, yy, zz = (pi*s/2 for s in [x,y,z])
        return Matrix([a*sin(2*xx),
                       b*cos(yy),
                       c*cos(zz),
                       cos(xx)*cos(yy)*cos(zz)])

def implementation():
  return '''
#include <vhtimpl.h>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

typedef struct {
  dReal a,b,c,scale;
} VHTCase_Exact;

static dErr VHTCaseSetFromOptions_Exact(VHTCase scase) {
  VHTCase_Exact *exc = scase->data;
  dErr err;

  dFunctionBegin;
  err = PetscOptionsHead("VHTCase_Exact options");dCHK(err); {
    err = PetscOptionsReal("-exact_a","First scale parameter","",exc->a,&exc->a,NULL);dCHK(err);
    err = PetscOptionsReal("-exact_b","Second scale parameter","",exc->b,&exc->b,NULL);dCHK(err);
    err = PetscOptionsReal("-exact_c","Third scale parameter","",exc->c,&exc->c,NULL);dCHK(err);
    err = PetscOptionsReal("-exact_scale","Overall scale parameter","",exc->scale,&exc->scale,NULL);dCHK(err);
  } err = PetscOptionsTail();dCHK(err);
  dFunctionReturn(0);
}
static dErr VHTCaseDestroy_Exact(VHTCase scase) {
  dErr err;

  dFunctionBegin;
  err = dFree(scase->data);dCHK(err);
  dFunctionReturn(0);
}
'''

def registerall(list):
    def register(sol):
        return 'err = VHTCaseRegister("%(name)s",VHTCaseCreate_%(name)s);dCHK(err);' % dict(name=sol.name)
    return '''
dErr VHTCaseRegisterAll_Exact(void)
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
    solutions = [Exact0(), Exact1()]
    with open('vhtexact.c', 'w') as fimpl:
        fimpl.write(implementation())
        for sol in solutions:
            fimpl.write(sol.solution_code())
            fimpl.write(sol.forcing_code())
            fimpl.write(sol.create_code())
        fimpl.write(registerall(solutions))
    print('Wrote vhtexact.c')
