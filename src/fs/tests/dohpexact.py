from __future__ import division

import sympy
from sympy import Matrix, ccode

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
def symbols(vars):
    sym = sympy.symbols(vars)
    try:                        # make sure that result is a list even if it only contains one entry
        sym[0]
        return sym
    except:
        return [sym]
def symbol3(var):
    return symbols(['%s[%d]' % (var, i) for i in range(3)])
def symbol33(var):
    return [symbol3('%s[%d]' % (var,i)) for i in range(3)]
I = sympy.eye(3)

class Exact:
    def __init__(self, fieldspec, name=None, model=None, param=None):
        self.name = name if name is not None else self.__class__.__name__
        self._fieldspec = fieldspec
        self.nfields = sum([x[1] for x in fieldspec])
        self._param = dict()
        self._model = dict()
        if model: self.model_add(model)
        if param: self.param_add(param)
        self.param_add('scale') # Always have small-displacement version available
    def fieldseek(self,i):
        for (name,count) in self._fieldspec:
            if i < count: return name, i
            i -= count
        raise KeyError
    def fieldname(self,i):
        return '%s[%d]' % self.fieldseek(i)
    def dfieldname(self,i):
        name, j = self.fieldseek(i//3)
        return 'd%s[%d]' % (name, j*3+i%3) # Index by flattened index
    def ffieldname(self,i):
        return 'f%s[%d]' % self.fieldseek(i)
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
        scale, = self.param_get('scale')
        return scale * self.solution(x[0], x[1], x[2], *args)
    def solution_gradient(self, x, *args):
        f = self.solution_scaled(x, *args)
        return grad(f,x)
    def residual(self, x, u, du):
        '''
        Converts from a weak form (scalar functional of the test
        functions) to residual form (coefficients of the test
        functions). The symbols v and dv disappear in differentiation.
        '''
        v = Matrix(symbols('v%d'%i for i in range(self.nfields)))
        dv = Matrix([symbols('v%d_%d'%(i,j) for j in range(3)) for i in range(self.nfields)])
        L = self.weak_homogeneous(x,u,du,v,dv)
        return testgrad(L,v), testgrad(L,dv)
