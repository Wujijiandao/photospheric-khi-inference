import numpy as np
from scipy.optimize import least_squares

def coefficients(kappa, epsilon):
    q=1-np.exp(-2*kappa)
    A=kappa*kappa*((epsilon*epsilon/4)*q-1)
    B=epsilon*kappa*q
    C=q+kappa*kappa-2*kappa-epsilon*epsilon*kappa*kappa*q/4
    return A,B,C

def mode(kappa, epsilon):
    A,B,C=coefficients(kappa,epsilon)
    D=B*B-4*A*C
    R=-B/(2*A)
    I=np.sqrt(max(0.0,-D))/(2*abs(A))
    return R,I,D

def invariants(d,epsilon,k1,k2):
    R1,I1,D1=mode(k1*d,epsilon)
    R2,I2,D2=mode(k2*d,epsilon)
    if I1<=0 or I2<=0 or D1>=0 or D2>=0:
        return np.array([np.nan,np.nan])
    return np.array([np.log(I1/I2),(R1-R2)/I1])

def invert_two_image_frame_modes(wavelength_km,growth_s_inv,vapp_km_s,
                                 d_bounds=(2.0,30.0),epsilon_bounds=(0.01,0.99)):
    """Recover d, epsilon, U0 and the common image-frame offset vc.
    Phase speeds must be signed along one common interface coordinate.
    """
    lam=np.asarray(wavelength_km,float)
    gamma=np.asarray(growth_s_inv,float)
    vapp=np.asarray(vapp_km_s,float)
    if len(lam)!=2:
        raise ValueError("Exactly two modes are required.")
    k=2*np.pi/lam
    g=gamma/k
    obs=np.array([np.log(g[0]/g[1]),(vapp[0]-vapp[1])/g[0]])

    def residual(x):
        th=invariants(x[0],x[1],k[0],k[1])
        if np.any(~np.isfinite(th)):
            return np.array([100.0,100.0])
        return th-obs

    fit=least_squares(residual,[10.0,0.5],
                      bounds=([d_bounds[0],epsilon_bounds[0]],
                              [d_bounds[1],epsilon_bounds[1]]))
    d,eps=fit.x
    RI=[mode(ki*d,eps)[:2] for ki in k]
    U0=np.mean([gi/(ki*I) for gi,ki,(R,I) in zip(gamma,k,RI)])
    vc=np.mean([vi-U0*R for vi,(R,I) in zip(vapp,RI)])
    return dict(d_km=float(d),epsilon=float(eps),U0_km_s=float(U0),
                vc_km_s=float(vc),residual_norm=float(np.linalg.norm(fit.fun)))
