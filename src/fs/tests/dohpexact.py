from __future__ import division, with_statement

import sympy
from sympy import Matrix, ccode
import itertools
#from IPython.Debugger import Tracer; itrace = Tracer()

def concat(llist):
    return itertools.chain(*llist)
def symdiff(f):
    "Differentiate a function of one variable, return a new function"
    X = sympy.Symbol('Xprivate')
    fX = f(X)                    # Symbolic expression for the function
    dfX = sympy.diff(fX,X)       # Symbolic expression for the derivative
    return sympy.lambdify(X,dfX) # Function that evaluates df
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
def asvector(A):
    m,n = A.shape
    return A.reshape(m*n,1)
def chaindiv(F, U, dUdX):
    assert(F.cols == dUdX.cols)
    assert(U.rows == dUdX.rows)
    assert(U.cols == 1)
    div = sympy.zeros((F.rows,1))
    for col in range(F.cols):
        div += grad(F[:,col], U) * dUdX[:,col]
    return div
def testgrad(A, v):
    return Matrix([[sympy.diff(A,v[i,j]) for j in range(v.shape[1])] for i in range(v.shape[0])])
def symbols(vars):
    sym = sympy.symbols(list(vars))
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
    def joinname(self, prefix, array, pair):
        return (prefix+('%s[%d]' if array else '%s_%d')) % pair
    def joinname2(self, prefix, array, trio):
        name, i, j = trio
        if array:
            return (prefix+'%s[%d]') % (name, i*3+j)
        else:
            return (prefix+'%s_%d_%d') % (name, i, j)
    def fieldname(self,i,prefix='',array=True):
        return self.joinname(prefix,array,self.fieldseek(i))
    def dfieldname(self,i,prefix='',array=True):
        name, j = self.fieldseek(i//3)
        return self.joinname2(prefix+'d',array,(name,j,i%3)) # Index by flattened index
    def ffieldname(self,i,prefix='',array=True):
        return self.joinname(prefix+'f',array,self.fieldseek(i))
    def fieldmatrices(self,prefix='',array=False):
        U = Matrix([self.fieldname(i,prefix,array) for i in range(self.nfields)])
        dU = Matrix([self.dfieldname(i*3+j,prefix,array) for i in range(self.nfields) for j in range(3)]).reshape(self.nfields,3)
        def d2name(dsym,k):
            'convert prefix_dsymname_i_j to prefix_d2symname_i_j_k'
            return prefix + 'd2' + dsym.name[len(prefix)+1:] + '_%d' % k
        d2U = Matrix([d2name(dsym,j) for dsym in asvector(dU) for j in range(3)]).reshape(self.nfields*3,3)
        return U, dU, d2U
    def _meta_add(self,attr,symnames):
        if isinstance(symnames,str): symnames = symnames.split()
        getattr(self,'_'+attr).update([(name,sympy.Symbol(name)) for name in symnames])
    def _meta_get(self,attr,symnames):
        if isinstance(symnames,str): symnames = symnames.split()
        return [getattr(self,'_'+attr)[name] for name in symnames]
    def model_add(self,symnames): return self._meta_add('model',symnames)
    def model_get(self,symnames): return self._meta_get('model',symnames)
    def param_add(self,symnames): return self._meta_add('param',symnames)
    def param_get(self,symnames): return self._meta_get('param',symnames)
    def param_decl(self):
        return 'const dReal ' + ','.join('%s = ctx->%s' % (v,v) for v in self._param.keys()) + ';'
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
    def residual_eval(self, x, u, du):
        "returns a list of C statements to evaluate the strong form of the residual"
        U, dU, d2U = self.fieldmatrices(prefix='X_',array=False)
        decls = map(vardeclaration,(U,dU,d2U))
        v, dv = self.residual(x,U,dU)
        preamble = zip(U,u) + zip(dU,du) + zip(d2U,grad(asvector(du),x))
        resid = v - chaindiv(dv,U,dU) - chaindiv(dv,asvector(dU),d2U)
        return decls, preamble, resid
    def residual_code(self, x, u, du):
        decls, preamble, resid = self.residual_eval(x, u, du)
        fnames = map(self.ffieldname,range(self.nfields))
        def mkstatement(pair):
            return ccode(pair[1], assign_to=pair[0])
        return decls, map(mkstatement,preamble) + map(mkstatement,zip(fnames,resid))

def vardeclaration(vars):
    lastname = ''
    def varname(v):
        head, sep, tail = v.name.partition('[')
        return head, bool(sep)
    def declname((key, iter)):
        head, sep = key
        n = len(list(iter))
        return '%s[%d]' % (head, sep) if sep else head
    gvars = itertools.groupby(vars, varname)
    decllist = map(declname,gvars)
    return 'dScalar dUNUSED ' + ','.join(decllist) + ';' # Set dUNUSED because we don't know if other expressions will use it

def utime(call,args=(),kwargs=dict()):
    import time
    start = time.time()
    result = call(*args,**kwargs)
    print call.__name__, time.time() - start
    return result
def uprofile(calls):
    import hotshot,hotshot.stats,tempfile
    logfile = tempfile.mktemp('.stats')
    try:
        prof = hotshot.Profile(logfile)
        for call,args,kwargs in calls:
            prof.runcall(call,args,kwargs)
        prof.close()
        stats = hotshot.stats.load(logfile)
        stats.print_stats(20)
    finally:
        os.remove(logfile)
