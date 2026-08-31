import numpy as np
from scipy.linalg import eig

def cheb(N):
    x=np.cos(np.pi*np.arange(N+1)/N)
    c=np.hstack([2.,np.ones(N-1),2.])*((-1.)**np.arange(N+1))
    X=np.tile(x,(N+1,1)).T
    dX=X-X.T
    D=np.outer(c,1/c)/(dX+np.eye(N+1))
    D=D-np.diag(np.sum(D,axis=1))
    return D,x

def solve(ka, epsilon, N=220, Lbox=12.0, density_width_ratio=1.0):
    """
    Generalized Rayleigh sensitivity model:
    U/U0=tanh(x/a), rho/rho0=1-epsilon*tanh(x/a_rho).
    x is nondimensionalized by a.
    """
    D,xi=cheb(N)
    x=Lbox*xi
    D1=D/Lbox; D2=D1@D1
    U=np.tanh(x)
    Up=1/np.cosh(x)**2
    Upp=-2*np.tanh(x)/(np.cosh(x)**2)
    ar=density_width_ratio
    rho=1-epsilon*np.tanh(x/ar)
    rhop=-(epsilon/ar)/(np.cosh(x/ar)**2)
    rlog=rhop/rho
    L=D2+np.diag(rlog)@D1-ka**2*np.eye(N+1)
    Q=Upp+rlog*Up
    A=np.diag(U)@L-np.diag(Q)
    B=L
    ii=np.arange(1,N)
    vals=eig(A[np.ix_(ii,ii)],B[np.ix_(ii,ii)],
             check_finite=False)[0]
    vals=vals[np.isfinite(vals)]
    c=vals[np.argmax(vals.imag)]
    return c, ka*c.imag
