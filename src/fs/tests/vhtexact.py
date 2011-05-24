#!/usr/bin/env python2
from __future__ import division

import sympy
from sympy import Matrix, ccode
from dohpexact import *

def SecondInvariant(Du):
    return 0.5*Du.dot(Du)
def const(x):
    def f(y): return x
    return f
def splice(a, b, x0, width, x, dx):
    from sympy import tanh
    ax = a(x)
    bx = b(x)
    dax = symdiff(a)(x)
    dbx = symdiff(b)(x)
    f = ax + (bx-ax) * (1+tanh((x-x0)/width)) / 2
    df = (dax
          + (dbx-dax) * (1+tanh((x-x0)/width)) / 2
          + (bx - ax) * (1-tanh((x-x0)/width)**2) / (2*width)) * dx
    return f, df

class VHTExact(Exact):
    def __init__(self, name=None, model=None, param='a b c'):
        if model is None:
            model = 'B0 Bomega R Q V T0 eps gamma0 pe beta_CC rhoi rhow T3 c_i Latent splice_delta k_T kappa_w'
        Exact.__init__(self, name=name, model=model, param=param, fieldspec=[('rhou',3), ('p',1), ('E',1)])
    def unpack(self, U, dU):
        rhou, p, E = U[:3,:], U[3], U[4]
        drhou, dp, dE = dU[:3,:], dU[3,:].T, dU[4,:].T
        return (rhou, p, E), (drhou, dp, dE)
    def solve_eqstate(self, rhou, p, E, drhou, dp, dE):
        "Uses the equation of state to solve for observable quantities and their derivatives"
        spdel, rhoi, rhow, T3, T0, c_i, beta_CC, L = self.model_get('splice_delta rhoi rhow T3 T0 c_i beta_CC Latent')
        T_m = T3 - beta_CC * p # melting temperature at current pressure
        e_m = c_i * (T_m - T0) # melting energy at current pressure
        # At this point, we should solve an implicit system for (e,rho,omega).  Unfortunately, doing so would require
        # solving an inhomogeneous system involving splice(), which I don't know how to do symbolically. Putting a
        # numeric rootfinder into this symbolic description would make manufactured solutions much more complex.
        # Therefore, we cheat by simply using ice density to remove kinetic energy and convert energy/volume to
        # energy/mass.
        rhotmp = rhoi          # Cheat
        e = (E - 1/(2*rhotmp) * rhou.dot(rhou)) / rhotmp
        de = (dE - 1/(rhotmp) * (rhou.T * drhou).T) / rhotmp
        def temp(e): return T0 + e/c_i # temperature in cold regime
        T, dT = splice(temp,const(T_m),e_m,spdel,e,de)
        def melt(e): return (e-e_m)/L # temperature in warm regime
        omega,domega  = splice(const(0),melt,e_m,spdel,e,de)
        rho = (1-omega)*rhoi + omega*rhow
        drho = (rhow-rhoi) * domega
        return (e,T,omega,rho), (de,dT,domega,drho)
    def eta(self, p, T, omega, gamma):       # Power law with Arrhenius relation
        from sympy import exp
        B0, Bomega, R, Q, V, T0, eps, gamma0, pe, beta_CC = self.model_get('B0 Bomega R Q V T0 eps gamma0 pe beta_CC')
        n = 1/(pe - 1)
        Tstar = T - beta_CC * p
        B = B0 * exp((Q*(T0 - Tstar) - p*V*T0) / (n*R*T0*Tstar)) * (1 + Bomega * omega)**(-1/n)
        return B0
        #return B * (eps**2 + gamma/gamma0)**((pe-2)/2)
    def weak_homogeneous(self, x, U, dU, V, dV):
        (rhou,p,E), (drhou,dp,dE) = self.unpack(U,dU)
        (e,T,omega,rho), (de,dT,domega,drho) = self.solve_eqstate(rhou,p,E,drhou,dp,dE)
        k_T, kappa_w, L = self.model_get('k_T kappa_w Latent')
        u = rhou / rho                    # total velocity
        du = (1/rho) * drhou - u * drho.T # total velocity gradient
        wmom = -kappa_w * domega          # momentum of water part in reference frame of ice, equal to mass flux
        ui = u - wmom/rho                 # ice velocity
        dui = du                          # We cheat again here. The second term should also be differentiated, but that would require second derivatives
        Dui = sym(dui)
        gamma = SecondInvariant(Dui)
        eta = self.eta(p, T, omega, gamma)
        heatflux = -k_T * dT + L * wmom
        (rhou_,p_,E_), (drhou_,dp_,dE_) = self.unpack(V,dV)
        conserve_momentum = -drhou_.dot(rhou*u.T - eta*Dui + p*I)
        conserve_mass     = -p_ * drhou.trace()
        conserve_energy   = -dE_.dot(E*ui + heatflux) - Dui.dot(eta*Dui)
        return conserve_momentum + conserve_mass + conserve_energy
    def create_prototype(self, definition=False):
        return 'dErr VHTCaseCreate_%(name)s(VHTCase case)%(term)s' % dict(name=self.name, term=('' if definition else ';'))
    def solution_prototype(self):
        return 'static dErr VHTCaseSolution_%(name)s(VHTCase scase,const dReal x[3],dScalar rhou[3],dScalar drhou[9],dScalar p[1],dScalar dp[3],dScalar E[1],dScalar dE[3])' % dict(name=self.name)
    def forcing_prototype(self):
        return 'static dErr VHTCaseForcing_%(name)s(VHTCase scase,const dReal x[3],dScalar frhou[3],dScalar fp[1],dScalar fE[1])' % dict(name=self.name)
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
            # decl, statements = self.residual_code(x,U,dU)
            # return decl + statements
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
    solutions = [Exact0()]
    with open('vhtexact.c', 'w') as fimpl:
        fimpl.write(implementation())
        for sol in solutions:
            fimpl.write(sol.solution_code())
            fimpl.write(sol.forcing_code())
            fimpl.write(sol.create_code())
        fimpl.write(registerall(solutions))
    print('Wrote vhtexact.c')
