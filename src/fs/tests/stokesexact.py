#!/usr/bin/env python2
from __future__ import division

import sympy
from sympy import Matrix, ccode
from dohpexact import *

def SecondInvariant(Du):
    return 0.5*Du.dot(Du)

class StokesExact(Exact):
    def __init__(self, name=None, model='A eps pe', param='a b c'):
        Exact.__init__(self, name=name, model=model, param=param, fieldspec=[('u',3), ('p',1)])
    def eta(self, gamma):       # Power law
        A, eps, pe = self.model_get('A eps pe')
        return A * (0.5*eps**2 + gamma)**((pe-2)/2)
    def weak_homogeneous(self, x, U, dU, V, dV):
        u, p = U[:3], U[3]
        Du = sym(dU[:3,:])
        v, q = V[:3], V[3]
        Dv = sym(dV[:3,:])
        gamma = SecondInvariant(Du)
        eta = self.eta(gamma)
        return eta*Dv.dot(Du) - q*Du.trace() - p*Dv.trace()
    def exact_prototype(self):
        return 'void %(name)s_Solution(const struct StokesExactCtx *ctx,const dReal x[3],dScalar u[3],dScalar du[9],dScalar p[1],dScalar dp[3])' % dict(name=self.name)
    def forcing_prototype(self):
        return 'void %(name)s_Forcing(const struct StokesExactCtx *ctx,const struct StokesRheology *rheo,const dReal x[3],dScalar fu[3],dScalar fp[1])' % dict(name=self.name)
    def exact_code(self):
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
  const dReal a = ctx->a,b = ctx->b,c = ctx->c,scale = ctx->scale;
  %(body)s
}''' % dict(prototype=self.exact_prototype(), body='\n  '.join(body()))
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
  const dReal a = ctx->a,b = ctx->b,c = ctx->c,scale = ctx->scale;
  const dReal A = rheo->A,eps = rheo->eps,pe = rheo->p;
  %(body)s
}
''' % dict(prototype=self.forcing_prototype(), body='\n  '.join(body()))

class StokesExact_0(StokesExact):
    'Classical 2D converging flow'
    def solution(self, x,y,z, a,b,c):
        from sympy import sin, cos, pi
        return Matrix([ a * sin(pi*x/2) * cos(pi*y/2),
                       -b * cos(pi*x/2) * sin(pi*y/2),
                        1 * (c-1) * z,
                        0.25*(cos(pi*x) + cos(pi*y)) + 10*(x+y)])
class StokesExact_1(StokesExact):
    'From Larin & Reusken, 2009, with pressure shifted by a constant'
    def solution(self, x,y,z, a,b,c):
        from sympy import sin,cos,pi
        xx, yy, zz = (pi*s/2 for s in [x,y,z])
        return Matrix([+a/3 * sin(xx) * sin(yy) * sin(zz),
                       -b/3 * cos(xx) * cos(yy) * sin(zz),
                       -c*2/3 * cos(xx) * sin(yy) * cos(zz),
                        1 + cos(xx) * sin(yy) * sin(zz)])
class StokesExact_2(StokesExact):
    def solution(self, x,y,z, a,b,c):
        return Matrix([a*z**3,
                       b*x**3,
                       c*y**3,
                       (1-x**2)*(1-y**2)*(1-z**2)])
class StokesExact_3(StokesExact):
    def solution(self, x,y,z, a,b,c):
        from sympy import sin,cos,pi
        xx, yy, zz = (pi*s/2 for s in [x,y,z])
        return Matrix([a*sin(2*xx),
                       b*cos(yy),
                       c*cos(zz),
                       cos(xx)*cos(yy)*cos(zz)])

if __name__ == "__main__":
    import pdb
    solutions = [StokesExact_0(), StokesExact_1(), StokesExact_2(), StokesExact_3()]
    with open('stokesexact.h', 'w') as fheader, open('stokesexact.c', 'w') as fimpl:
        fheader.write('#include <dohptype.h>\n\n')
        fheader.write('struct StokesRheology {dReal A,eps,p;};\n')
        fheader.write('struct StokesExactCtx {dReal a,b,c,scale;};\n')
        fimpl.write('\n'.join(['#include "stokesexact.h"', '#include <math.h>', '#ifndef M_PI', '#  define M_PI 3.14159265358979323846', '#endif', '\n']))
        for sol in solutions:
            fheader.write(sol.exact_prototype() + ';\n')
            fheader.write(sol.forcing_prototype() + ';\n')
            fimpl.write(sol.exact_code())
            fimpl.write(sol.forcing_code())
    print('Wrote stokesexact.h and stokesexact.c')
