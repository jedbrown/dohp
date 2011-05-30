#!/usr/bin/env python2
from __future__ import division, with_statement

import sympy
from sympy import Matrix, Symbol, ccode, S
from dohpexact import *

# Some coupling terms make the symbolic expressions explode in size, which takes too long and crashes the compiler. Such
# terms are multiplied by MASK so that they can be easily activated by setting MASK=1. There are a couple other places,
# marked by comments containing the word "cheat", where the fully consistent constitutive relation requires solving a
# numerical rootfinding problem. Since we cannot symbolically differentiate an implicit solve, we are shorting it out in
# a physically plausible way, but one that is inconsistent for large amplitude moisture.
MASK = 0
MASK2 = 0

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
def normal(x,w):
    "normal distribution with standard deviation w"
    return exp(-x**2/(2*w**2))/sqrt(2*pi*w**2)
def separate(x,w):
    """Splits the function f(x)=x into a negative and positive part. The derivatives on each side are the error function
    and a Gaussian, so unlike most expressions, these get nicer when they are differentiated.
    neg = integrate(0.5+integrate(normal(x,a),x),x)
    pos = integrate(0.5-integrate(normal(x,a),x),x)
    """
    from sympy import erf,exp,pi,sqrt
    # This explicit CSE does not seem to make a performance difference with sympy, but it probably would if we use AD
    theerf = erf(x/sqrt(2*w**2))
    theexp = exp(-x**2/(2*w**2))
    neg = 0.5*x*(1 - theerf) - w*theexp/sqrt(2*pi)
    pos = 0.5*x*(1 + theerf) + w*theexp/sqrt(2*pi)
    neg1 = 0.5*(1 - theerf) # first derivwtives
    pos1 = 0.5*(1 + theerf)
    return neg,pos,neg1,pos1

class VHTExact(Exact):
    def __init__(self, name=None, model=None, param='a b c d e'):
        if model is None:
            model = ' '.join(['B0 Bomega R Q V T0 p0 eps gamma0 pe beta_CC rhoi rhow T3 c_i Latent splice_delta k_T kappa_w gravity',
                              'Kstab mask_momtrans mask_rho mask_Ep',
                              ])
        Exact.__init__(self, name=name, model=model, param=param, fieldspec=[('rhou',3), ('p',1), ('E',1)])
    def unpack(self, U, dU):
        rhou, p, E = U[:3,:], U[3], U[4]
        drhou, dp, dE = dU[:3,:], dU[3,:].T, dU[4,:].T
        return (rhou, p, E), (drhou, dp, dE)
    def solve_eqstate(self, rhou, p, E, drhou, dp, dE):
        "Uses the equation of state to solve for observable quantities and their derivatives"
        spdel, rhoi, rhow, T3, T0, p0, c_i, beta_CC, L = self.model_get('splice_delta rhoi rhow T3 T0 p0 c_i beta_CC Latent')
        T_m = T3 - beta_CC * (p+p0) # melting temperature at current pressure
        e_m = c_i * (T_m - T0) # melting energy at current pressure
        # At this point, we should solve an implicit system for (e,rho,omega).  Unfortunately, doing so would require
        # solving an inhomogeneous system involving separate(), is problematic to do symbolically. Putting a
        # numeric rootfinder into this symbolic description would make manufactured solutions much more complex.
        # Therefore, we cheat by simply using ice density to remove kinetic energy and convert energy/volume to
        # energy/mass. Note that with algorithmic differentiation, the implicit solve would not be a problem.
        rhotmp = rhoi          # Cheat
        e  = (E - MASK*1/(2*rhotmp) * rhou.dot(rhou)) / rhotmp
        de = (dE - MASK*1/(rhotmp) * (rhou.T * drhou).T) / rhotmp
        G   = e - e_m
        G1p = c_i * beta_CC
        G1E = 1/rhotmp;
        e4T, e4omega, e4T1, e4omega1 = separate(G, spdel) # Decompose the energy G into a thermal and moisture part, with derivatives
        T      = T0 + (e_m+e4T)/c_i
        dT     = -beta_CC*dp + e4T1*(G1p*dp + G1E*dE)/c_i
        def meltfraction(e):
            'Convert energy/mass to melt fraction'
            x = Symbol('x')
            #omega = rhoi*x / (rhow*L - (rhow-rhoi)*x)
            omega = rhoi/(rhow*L)*x # Cheat: should use expression above
            omega1 = sympy.diff(omega,x)
            return omega.subs(x,e), omega1.subs(x,e)
        omega, omega1 = meltfraction(e4omega1)
        domega = omega1 * e4omega1 * (G1p*dp + G1E*dE)
        if MASK2 == 0: # Shameful, but SymPy chokes
            rho = rhoi
            drho = 0*domega
        else:
            m = self.model_get('mask_rho')
            rho = (1-m*omega)*rhoi + m*omega*rhow
            drho = (rhow-rhoi) * m*domega
        return (e,T,omega,rho), (de,dT,domega,drho)
    def eta(self, p, T, omega, gamma):       # Power law with Arrhenius relation
        from sympy import exp
        B0, Bomega, R, Q, V, T0, p0, eps, gamma0, pe, beta_CC = self.model_get('B0 Bomega R Q V T0 p0 eps gamma0 pe beta_CC')
        n = 1/(pe - 1)
        Tstar = T - beta_CC * (p+p0)
        B = B0 * exp((Q - (p+p0)*V) / (n*R*Tstar) - Q/(n*R*T0)) * (1 + Bomega * omega)**(-1/n)
        return B * (eps**2 + gamma/gamma0)**((pe-2)/2)
    def weak_homogeneous(self, x, U, dU, V, dV):
        (rhou,p,E), (drhou,dp,dE) = self.unpack(U,dU)
        (e,T,omega,rho), (de,dT,domega,drho) = self.solve_eqstate(rhou,p,E,drhou,dp,dE)
        k_T, kappa_w, Kstab, L, grav, rhoi, p0 = self.model_get('k_T kappa_w Kstab Latent gravity rhoi p0')
        mask_momtrans, mask_Ep = self.model_get('mask_momtrans mask_Ep')
        gravvec = grav*Matrix([0,0,1])
        u = rhou / rho                    # total velocity
        du = (1/rho) * drhou - MASK*u * drho.T # total velocity gradient
        dui = du                          # We cheat again here. The second term should also be differentiated, but that would require second derivatives
        Dui = sym(dui)
        gamma = SecondInvariant(Dui)
        eta = self.eta(p, T, omega, gamma)
        heatflux = -k_T * dT - L * (1-omega)*rhoi / rho * kappa_w * domega;
        (rhou_,p_,E_), (drhou_,dp_,dE_) = self.unpack(V,dV)
        conserve_momentum = -drhou_.dot(mask_momtrans*rhou*u.T - eta*Dui + p*I) - rhou_.dot(rho*gravvec)
        conserve_mass     = -p_ * drhou.trace()
        conserve_energy   = -dE_.dot((E+mask_Ep*(p+p0))*u + heatflux - Kstab*dE) -E_*(eta*Dui.dot(Dui) + rhou.dot(gravvec))
        return conserve_momentum + conserve_mass + conserve_energy
    def create_prototype(self, definition=False):
        return 'dErr VHTCaseCreate_%(name)s(VHTCase case)%(term)s' % dict(name=self.name, term=('' if definition else ';'))
    def solution_prototype(self):
        return 'static dErr VHTCaseSolution_%(name)s(VHTCase scase,const dReal x[3],dScalar rhou[3],dScalar drhou[9],dScalar p[1],dScalar dp[3],dScalar E[1],dScalar dE[3])' % dict(name=self.name)
    def forcing_prototype(self):
        return 'static dErr VHTCaseForcing_%(name)s(VHTCase scase,const dReal x[3],dScalar frhou[3],dScalar fp[1],dScalar fE[1])' % dict(name=self.name)
    def solution_code(self):
        from sympy.abc import a,b,c,d,e
        x = Matrix(symbol3('x'))
        def body():
            for (i,expr) in enumerate(self.solution_scaled(x,a,b,c,d,e)):
                yield ccode(expr, assign_to=self.fieldname(i))
            for (i,expr) in enumerate(self.solution_gradient(x,a,b,c,d,e)):
                yield ccode(expr, assign_to=self.dfieldname(i))
        return '''
%(prototype)s
{
  const VHTCase_Exact *ctx = scase->data;
  %(param)s
  %(body)s
  return 0;
}''' % dict(prototype=self.solution_prototype(), param=self.param_decl(), body='\n  '.join(body()))
    def forcing_code(self):
        from sympy.abc import a,b,c,d,e
        x = Matrix(symbol3('x'))
        def body():
            U = self.solution_scaled(x,a,b,c,d,e)
            dU = self.solution_gradient(x,a,b,c,d,e)
            if True:
                decl, statements = self.residual_code(x,U,dU)
                return decl + statements
            else:
                V,dV = self.residual(x,U,dU)
                return [ccode(V[i] - divergence(row,x), assign_to=self.ffieldname(i)) for i,row in enumerate(rows(dV))]
        return '''
%(prototype)s
{
  const VHTCase_Exact *ctx = scase->data;
  %(param)s
  const dUNUSED dReal %(rheo)s;
  %(body)s
  return 0;
}
''' % dict(prototype=self.forcing_prototype(),
           param=self.param_decl(),
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
    def solution(self, x,y,z, a,b,c,d,e):
        from sympy import sin, cos, pi
        return Matrix([ a * sin(pi*x/2) * cos(pi*y/2),
                       -b * cos(pi*x/2) * sin(pi*y/2),
                        1 * (c-1) * z,
                        d*0.25*(cos(pi*x) + cos(pi*y)) + 2*(x+y),
                        e*sin(pi*(x+y+z*z)/2) * cos(pi*(x*x+y+z)/2)])
class Exact1(VHTExact):
    'From Larin & Reusken, 2009, with pressure shifted by a constant'
    def solution(self, x,y,z, a,b,c,d,e):
        from sympy import sin,cos,pi,tanh,exp
        xx, yy, zz = (pi*s/2 for s in [x,y,z])
        r2 = xx**2 + yy**2 + zz**2
        return Matrix([+a/3 * sin(xx) * sin(yy) * sin(zz),
                       -b/3 * cos(xx) * cos(yy) * sin(zz),
                       -c*2/3 * cos(xx) * sin(yy) * cos(zz),
                        1 + d * cos(xx) * sin(yy) * sin(zz),
                        e*sin(pi*(x+y+z*z)/2) * cos(pi*(x*x+y+z)/2)])
class Exact2(VHTExact):
    'From Larin & Reusken, 2009, pressure shifted by a constant, spherical energy contours with cold spot in the center'
    def solution(self, x,y,z, a,b,c,d,e):
        from sympy import sin,cos,pi,tanh,exp
        xx, yy, zz = (pi*s/2 for s in [x,y,z])
        r2 = xx**2 + yy**2 + zz**2
        return Matrix([+a/3 * sin(xx) * sin(yy) * sin(zz),
                       -b/3 * cos(xx) * cos(yy) * sin(zz),
                       -c*2/3 * cos(xx) * sin(yy) * cos(zz),
                        1 + d * cos(xx) * sin(yy) * sin(zz),
                        e * 10 * r2 * exp(-4*r2)])

def implementation():
  return '''
#include <vhtimpl.h>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

typedef struct {
  dReal a,b,c,d,e,scale;
} VHTCase_Exact;

static dErr VHTCaseSetFromOptions_Exact(VHTCase scase) {
  VHTCase_Exact *exc = scase->data;
  dErr err;

  dFunctionBegin;
  err = PetscOptionsHead("VHTCase_Exact options");dCHK(err); {
    err = PetscOptionsReal("-exact_a","First scale parameter","",exc->a,&exc->a,NULL);dCHK(err);
    err = PetscOptionsReal("-exact_b","Second scale parameter","",exc->b,&exc->b,NULL);dCHK(err);
    err = PetscOptionsReal("-exact_c","Third scale parameter","",exc->c,&exc->c,NULL);dCHK(err);
    err = PetscOptionsReal("-exact_d","Pressure scale parameter","",exc->d,&exc->d,NULL);dCHK(err);
    err = PetscOptionsReal("-exact_e","Energy scale parameter","",exc->e,&exc->e,NULL);dCHK(err);
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
