import numpy as np
from scipy.optimize import minimize_scalar, brentq

def coefficients(kappa, epsilon):
    q = 1.0 - np.exp(-2.0*kappa)
    A = kappa**2*((epsilon**2/4.0)*q - 1.0)
    B = epsilon*kappa*q
    C = q + kappa**2 - 2.0*kappa - epsilon**2*kappa**2*q/4.0
    return A,B,C

def mode(kappa, epsilon):
    A,B,C = coefficients(kappa, epsilon)
    D = B*B - 4*A*C
    R = -B/(2*A)
    I = np.sqrt(max(0.0,-D))/(2*abs(A))
    return R,I,D

def fastest_mode(epsilon):
    r = minimize_scalar(lambda k: -k*mode(k,epsilon)[1],
                        bounds=(1e-6,1.8),method="bounded")
    k = r.x
    R,I,_ = mode(k,epsilon)
    G = k*I
    return {
        "epsilon":epsilon, "kappa_star":k, "lambda_over_d":2*np.pi/k,
        "gamma_d_over_U0":G, "vph_over_U0":R,
        "chi":R/(2*np.pi*I) if I else 0.0
    }
