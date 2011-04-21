#!/usr/bin/env python2
from __future__ import division

import sympy
from sympy import Matrix, symbols, ccode
from dohpexact import *

class ElastExact(Exact):
    def __init__(self, name=None, model='lambda mu', param='a b c'):
        Exact.__init__(self, name=name, model=model, param=param, nfields=3)
    def PiolaKirchoff2(self, E):
        lmbda, mu = self.model_get('lambda mu')
        return lmbda*E.trace()*I + 2*mu*E # St. Venant-Kirchoff model
    def weak_homogeneous(self, x, u, du, v, dv):
        H = du                          # Displacement gradient
        F = I + H                       # Deformation gradient
        C = F.T * F                     # Cauchy-Green tensor
        E = 0.5 * (C - I)               # Green-Lagrange tensor
        E_alt = 0.5 * (H + H.T + H.T*H) # Alternative expression
        S = self.PiolaKirchoff2(E)      # Second Piola-Kirchoff tensor
        Pi = F * S                      # First Piola-Kirchoff tensor
        #L0 = 0.5*(dv + dv.T + dv.T*H + H.T*dv).dot(S) # Obviously symmetric form
        #L1 = (F.T * dv).dot(S)                        # Equivalent because S is symmetric
        L2 = dv.dot(Pi)                               # Move F to the other side of the contraction by identity (obvious in index notation)
        return L2
    def exact_prototype(self):
        return 'void %(name)s_Solution(const struct ElastExactCtx *ctx,const dReal x[3],dScalar u[3],dScalar du_flat[9])' % dict(name=self.name)
    def forcing_prototype(self):
        return 'void %(name)s_Forcing(const struct ElastExactCtx *ctx,const struct ElastParam *prm,const dReal x[3],dScalar f[3])' % dict(name=self.name)
    def exact_code(self):
        from sympy.abc import a,b,c
        x = Matrix(symbol3('x'))
        def body():
            for (i,expr) in enumerate(self.solution_scaled(x,a,b,c)):
                yield ccode(expr, assign_to='u[%d]'%i)
            for (i,expr) in enumerate(self.solution_gradient(x,a,b,c)):
                yield ccode(expr, assign_to='du_flat[%d]'%i)
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
            u = self.solution_scaled(x,a,b,c)
            du = self.solution_gradient(x,a,b,c)
            v,dv = self.residual(x,u,du)
            for i,row in enumerate(rows(dv)):
                yield ccode(v[i] - divergence(row,x), assign_to='f[%d]'%i)
        return '''
%(prototype)s
{
  const dReal a = ctx->a,b = ctx->b,c = ctx->c,scale = ctx->scale;
  const dReal lambda = prm->lambda,mu = prm->mu;
  %(body)s
}
''' % dict(prototype=self.forcing_prototype(), body='\n  '.join(body()))

class ElastExact_0(ElastExact):
    'Product of transcendental functions'
    def solution(self, x,y,z, a,b,c):
        from sympy import cos,sin,exp,cosh,sinh,tanh,log
        return Matrix([cos(a*x) * exp(b*y) * z + sin(c*z),
                       sin(a*x) * tanh(b*y) + x * cosh(c*z),
                       exp(a*x) * sinh(b*y) + y * log(1+(c*z)**2)])
class ElastExact_1(ElastExact):
    '''
    Displacement-only manufactured solution used in chamberland2010comparison
    title={Comparison of the performance of some finite element discretizations for large deformation elasticity problems}
    author={Chamberland, {\'E}. and Fortin, A. and Fortin, M.}
    year={2010}
    The solution in the paper uses a=b=c=1
    '''
    def solution(self, x,y,z, a,b,c):
        return Matrix([x**4 + 2*a*y*z/5,
                       b*y**4 + 2*x*z/5,
                       z**4/10 - c*2*x*y*z/5])
class ElastExact_2(ElastExact):
    def solution(self, x,y,z, a,b,c):
        return Matrix([a*z**3, b*z**3, c*0.01*z])

if __name__ == "__main__":
    solutions = [ElastExact_0(), ElastExact_1(), ElastExact_2()]
    with open('elastexact.h', 'w') as fheader, open('elastexact.c', 'w') as fimpl:
        fheader.write('#include <dohptype.h>\n\n')
        fheader.write('struct ElastParam {dReal lambda,mu;};\n')
        fheader.write('struct ElastExactCtx {dReal a,b,c,scale;};\n')
        fimpl.write('#include "elastexact.h"\n')
        fimpl.write('#include <math.h>\n\n')
        for sol in solutions:
            fheader.write(sol.exact_prototype() + ';\n')
            fheader.write(sol.forcing_prototype() + ';\n')
            fimpl.write(sol.exact_code())
            fimpl.write(sol.forcing_code())
    print('Wrote elastexact.h and elastexact.c')
