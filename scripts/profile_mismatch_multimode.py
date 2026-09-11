#!/usr/bin/env python3
"""Three-mode profile-mismatch test for the image-frame KHI inverse.

Generates three KH modes from smooth variable-density tanh, erf, and arctan
profiles, fits the outer two modes with the finite-slab image-frame inverse,
and uses the middle mode as an internal closure test.
"""
import numpy as np
import pandas as pd
from scipy.linalg import eig
from scipy.optimize import brentq, least_squares
from scipy.special import erf, erfinv
from pathlib import Path

OUT = Path(__file__).resolve().parent / "profile_mismatch_multimode_test.csv"


def coefficients(kappa, epsilon):
    q = 1.0 - np.exp(-2.0*kappa)
    A = kappa**2*((epsilon**2/4.0)*q - 1.0)
    B = epsilon*kappa*q
    C = q + kappa**2 - 2.0*kappa - epsilon**2*kappa**2*q/4.0
    return A, B, C


def slab_mode(kappa, epsilon):
    A, B, C = coefficients(kappa, epsilon)
    D = B*B - 4*A*C
    R = -B/(2*A)
    I = np.sqrt(max(0.0, -D))/(2*abs(A))
    return R, I, D


def slab_cutoff(epsilon):
    xs = np.linspace(1e-4, 3.0, 3000)
    ds = np.array([slab_mode(x, epsilon)[2] for x in xs])
    ix = np.where((ds[:-1] < 0) & (ds[1:] >= 0))[0]
    j = ix[-1]
    return brentq(lambda x: slab_mode(x, epsilon)[2], xs[j], xs[j+1])


def cheb(N):
    x = np.cos(np.pi*np.arange(N+1)/N)
    c = np.hstack([2., np.ones(N-1), 2.])*((-1.)**np.arange(N+1))
    X = np.tile(x, (N+1, 1)).T
    dX = X-X.T
    D = np.outer(c, 1/c)/(dX+np.eye(N+1))
    D = D-np.diag(np.sum(D, axis=1))
    return D, x


def profile_shape(x, family):
    if family == "tanh":
        f = np.tanh(x)
        fp = 1/np.cosh(x)**2
        fpp = -2*np.tanh(x)/(np.cosh(x)**2)
        width = 2*np.arctanh(0.8)
    elif family == "erf":
        f = erf(x)
        fp = 2/np.sqrt(np.pi)*np.exp(-x*x)
        fpp = -4*x/np.sqrt(np.pi)*np.exp(-x*x)
        width = 2*erfinv(0.8)
    elif family == "arctan":
        f = 2/np.pi*np.arctan(x)
        fp = 2/np.pi/(1+x*x)
        fpp = -4*x/np.pi/(1+x*x)**2
        width = 2*np.tan(0.4*np.pi)
    else:
        raise ValueError(family)
    return f, fp, fpp, width


def smooth_c(ka, epsilon, family, N=100):
    Lbox = 20.0 if family == "arctan" else 12.0
    D, xi = cheb(N)
    x = Lbox*xi
    D1 = D/Lbox
    D2 = D1@D1
    f, fp, fpp, _ = profile_shape(x, family)
    rho = 1-epsilon*f
    rlog = (-epsilon*fp)/rho
    Lop = D2 + np.diag(rlog)@D1 - ka**2*np.eye(N+1)
    A = np.diag(f)@Lop - np.diag(fpp+rlog*fp)
    B = Lop
    ii = np.arange(1, N)
    vals = eig(A[np.ix_(ii, ii)], B[np.ix_(ii, ii)], check_finite=False)[0]
    vals = vals[np.isfinite(vals)]
    vals = vals[vals.imag > 1e-10]
    if len(vals) == 0:
        raise RuntimeError(f"No unstable eigenvalue for {family}, ka={ka}")
    return vals[np.argmax(vals.imag)]


def smooth_forward(family, fracs=(0.30, 0.60, 0.90),
                   d=10.0, U0=1.5, epsilon=0.6, vc=1.2):
    kc = slab_cutoff(epsilon)
    k = np.array(fracs)*kc/d
    lam = 2*np.pi/k
    _, _, _, width = profile_shape(np.array([0.0]), family)
    a = d/width
    gamma, vapp = [], []
    for ki in k:
        c = smooth_c(ki*a, epsilon, family)
        # Orient the along-interface coordinate so propagation has the same
        # sign convention as the slab branch used in the inverse.
        R = -c.real
        I = c.imag
        gamma.append(ki*U0*I)
        vapp.append(vc+U0*R)
    return lam, np.array(gamma), np.array(vapp)


def invariants_obs(lam, gamma, vapp):
    k = 2*np.pi/np.asarray(lam)
    g = np.asarray(gamma)/k
    return np.array(
        [np.log(g[0]/g[j]) for j in range(1, len(k))] +
        [(vapp[0]-vapp[j])/g[0] for j in range(1, len(k))]
    )


def invariants_slab_n(d, epsilon, k):
    RI = [slab_mode(ki*d, epsilon)[:2] for ki in k]
    if any(I <= 1e-10 for R, I in RI):
        return None
    R = np.array([z[0] for z in RI])
    I = np.array([z[1] for z in RI])
    return np.array(
        [np.log(I[0]/I[j]) for j in range(1, len(k))] +
        [(R[0]-R[j])/I[0] for j in range(1, len(k))]
    )


def fit_slab_modes(lam, gamma, vapp, indices):
    lam = np.asarray(lam)[indices]
    gamma = np.asarray(gamma)[indices]
    vapp = np.asarray(vapp)[indices]
    k = 2*np.pi/lam
    obs = invariants_obs(lam, gamma, vapp)

    def residual(x):
        th = invariants_slab_n(x[0], x[1], k)
        if th is None or np.any(~np.isfinite(th)):
            return np.ones_like(obs)*100
        return th-obs

    fit = least_squares(residual, [10.0, 0.6],
                        bounds=([2.0, 0.01], [30.0, 0.99]), max_nfev=300)
    d, epsilon = fit.x
    RI = np.array([slab_mode(ki*d, epsilon)[:2] for ki in k])
    I, R = RI[:, 1], RI[:, 0]
    g = gamma/k
    U0 = float(np.sum(I*g)/np.sum(I*I))
    vc = float(np.mean(vapp-U0*R))
    return d, epsilon, U0, vc, float(np.linalg.norm(fit.fun))


def main():
    rows = []
    for family in ["tanh", "erf", "arctan"]:
        lam, gamma, vapp = smooth_forward(family)
        d, e, U0, vc, _ = fit_slab_modes(lam, gamma, vapp, [0, 2])
        kmid = 2*np.pi/lam[1]
        Rm, Im, _ = slab_mode(kmid*d, e)
        gamma_pred = kmid*U0*Im
        vapp_pred = vc+U0*Rm
        d3, e3, U3, vc3, res3 = fit_slab_modes(lam, gamma, vapp, [0, 1, 2])
        rows.append({
            "profile": family,
            "truth_d_km": 10.0, "truth_epsilon": 0.6,
            "truth_U0_km_s": 1.5, "truth_vc_km_s": 1.2,
            "lambda1_km": lam[0], "lambda2_km": lam[1], "lambda3_km": lam[2],
            "outer_pair_d_km": d, "outer_pair_epsilon": e,
            "outer_pair_U0_km_s": U0, "outer_pair_vc_km_s": vc,
            "d_bias_percent": 100*(d/10.0-1),
            "epsilon_bias_percent": 100*(e/0.6-1),
            "U0_bias_percent": 100*(U0/1.5-1),
            "vc_bias_percent": 100*(vc/1.2-1),
            "middle_gamma_true_s-1": gamma[1],
            "middle_gamma_pred_s-1": gamma_pred,
            "middle_gamma_error_percent": 100*(gamma_pred/gamma[1]-1),
            "middle_vapp_true_km_s": vapp[1],
            "middle_vapp_pred_km_s": vapp_pred,
            "middle_vapp_error_km_s": vapp_pred-vapp[1],
            "three_mode_fit_d_km": d3, "three_mode_fit_epsilon": e3,
            "three_mode_fit_U0_km_s": U3, "three_mode_fit_vc_km_s": vc3,
            "three_mode_invariant_residual": res3,
        })
    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
