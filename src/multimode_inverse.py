import numpy as np
from scipy.optimize import least_squares


def coefficients(k,e):
    q=1-np.exp(-2*k)
    A=k*k*((e*e/4)*q-1)
    B=e*k*q
    C=q+k*k-2*k-e*e*k*k*q/4
    return A,B,C


def mode(k,e):
    A,B,C=coefficients(k,e)
    D=B*B-4*A*C
    R=-B/(2*A)
    I=np.sqrt(max(0.0,-D))/(2*abs(A))
    return R,I


def Z(k,e):
    R,I=mode(k,e)
    return R/I if I>1e-12 else np.nan


def invert_two_modes(wavelengths_km,growth_s_inv,vph_km_s,
                     d_bounds=(2.0,30.0),e_bounds=(0.02,0.98)):
    wavelengths_km=np.asarray(wavelengths_km,float)
    growth_s_inv=np.asarray(growth_s_inv,float)
    vph_km_s=np.asarray(vph_km_s,float)
    k=2*np.pi/wavelengths_km
    zeta=k*vph_km_s/growth_s_inv

    def residual(x):
        d,e=x
        zz=np.array([Z(ki*d,e) for ki in k])
        if np.any(~np.isfinite(zz)):
            return np.full(2,100.0)
        return zz-zeta

    fit=least_squares(residual,[10.0,0.5],
                      bounds=([d_bounds[0],e_bounds[0]],
                              [d_bounds[1],e_bounds[1]]))
    d,e=fit.x
    U=[]
    for ki,gi,vi in zip(k,growth_s_inv,vph_km_s):
        R,I=mode(ki*d,e)
        U.extend([vi/R,gi/(ki*I)])
    return {"d_km":d,"epsilon":e,"U0_km_s":float(np.median(U)),
            "residual_norm":float(np.linalg.norm(fit.fun))}
