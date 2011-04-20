#!/usr/bin/env python2
from __future__ import division

import sympy
from sympy import Matrix, symbols, ccode

def rows(A):
    return [A[i,:] for i in range(A.rows)]
def grad(U, X):
    return U.jacobian(X)
def symgrad(U, X):
    return sym(grad(U,X))
def sym(A):
    return 0.5*(A + A.T)
def divergence(F, X):
    assert(F.cols == X.rows)
    return sum(sympy.diff(F[i],X[i]) for i in range(X.rows))
def testgrad(A, v):
    return Matrix([[sympy.diff(A,v[i,j]) for j in range(v.shape[1])] for i in range(v.shape[0])])
def symbol3(var):
    return symbols(['%s[%d]' % (var, i) for i in range(3)])
def symbol33(var):
    return [symbol3('%s[%d]' % (var,i)) for i in range(3)]
I = sympy.eye(3)

class Elast:
    def __init__(self, name='ElastExact_X', param=None):
        self.name = name
        self._param = dict()
        self._model = dict()
        self.model_add('lambda mu scale')
        if param: self.param_add(param)
    def _meta_add(self,attr,symnames):
        if isinstance(symnames,str): symnames = symnames.split()
        getattr(self,'_'+attr).update([(name,sympy.symbols([name])) for name in symnames])
    def _meta_get(self,attr,symnames):
        if isinstance(symnames,str): symnames = symnames.split()
        return [getattr(self,'_'+attr)[name] for name in symnames]
    def model_add(self,symnames): return self._meta_add('model',symnames)
    def model_get(self,symnames): return self._meta_get('model',symnames)
    def param_add(self,symnames): return self._meta_add('param',symnames)
    def param_get(self,symnames): return self._meta_get('param',symnames)
    def solution(self, x, *args):
        raise NotImplementedError
    def solution_scaled(self, x, *args):
        scale, = self.model_get('scale')
        return scale * self.solution(x[0], x[1], x[2], *args)
    def solution_gradient(self, x, *args):
        f = self.solution_scaled(x, *args)
        return grad(f,x)
    def PiolaKirchoff2(self, E):
        lmbda, mu = self.model_get('lambda mu')
        return lmbda*E.trace()*I + 2*mu*E # Saint-Venant Kirchoff model
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
    def residual(self, x, u, du):
        '''
        Converts from a weak form (scalar functional of the test
        functions) to residual form (coefficients of the test
        functions).
        '''
        v  = Matrix(symbol3('v'))
        dv = Matrix(symbol33('dv'))
        L = self.weak_homogeneous(x,u,du,v,dv)
        return testgrad(L,v), testgrad(L,dv)
    def code(self):
        x  = Matrix(symbol3('x'))
        u  = Matrix(symbol3('u'))
        du = Matrix(symbol33('du'))
        v, dv = self.residual(x,u,du)
        for (i,x) in enumerate(dv):
            yield ccode(x, assign_to='Dv_flat[%d]'%i)
        # print '### Entries in Jacobian'
        # for part in sympy.cse(testgrad(test_dv[0,0],du)):
        #     for x in part:
        #         yield x
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

class Elast_0(Elast):
    'Product of transcendental functions'
    def __init__(self): Elast.__init__(self, name='ElastExact_0', param='a b c')
    def solution(self, x,y,z, a,b,c):
        from sympy import cos,sin,exp,cosh,sinh,tanh,log
        return Matrix([cos(a*x) * exp(b*y) * z + sin(c*z),
                       sin(a*x) * tanh(b*y) + x * cosh(c*z),
                       exp(a*x) * sinh(b*y) + y * log(1+(c*z)**2)])
class Elast_1(Elast):
    '''
    Displacement-only manufactured solution used in chamberland2010comparison
    title={Comparison of the performance of some finite element discretizations for large deformation elasticity problems}
    author={Chamberland, {\'E}. and Fortin, A. and Fortin, M.}
    year={2010}
    The solution in the paper uses a=b=c=1
    '''
    def __init__(self): Elast.__init__(self, name='ElastExact_1', param='a b c')
    def solution(self, x,y,z, a,b,c):
        return Matrix([x**4 + 2*a*y*z/5,
                       b*y**4 + 2*x*z/5,
                       z**4/10 - c*2*x*y*z/5])
class Elast_2(Elast):
    def __init__(self): Elast.__init__(self, name='ElastExact_2', param='a b c')
    def solution(self, x,y,z, a,b,c):
        return Matrix([a*z**3, b*z**3, c*0.01*z])

if __name__ == "__main__":
    solutions = [Elast_0(), Elast_1(), Elast_2()]
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
