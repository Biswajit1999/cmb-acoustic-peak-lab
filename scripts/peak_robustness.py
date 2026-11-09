"""Local CMB binned-spectrum peak diagnostics, not cosmological inference."""
import numpy as np

def quadratic_peak(ell, power, sigma, covariance=None):
    """Weighted quadratic vertex and linearized covariance uncertainty."""
    x, y, s = [np.asarray(z, float) for z in (ell,power,sigma)]
    if x.ndim != 1 or len(x) < 3 or x.shape != y.shape or x.shape != s.shape:
        raise ValueError("at least three aligned bins required")
    if not all(np.all(np.isfinite(z)) for z in (x,y,s)) or np.any(s<=0) or np.any(np.diff(x)<=0):
        raise ValueError("finite sorted multipoles and positive errors required")
    u=x-np.mean(x)
    X=np.column_stack([u*u,u,np.ones_like(u)])
    C=np.diag(s*s) if covariance is None else np.asarray(covariance,float)
    if C.shape != (len(x),len(x)) or not np.all(np.isfinite(C)) or not np.allclose(C,C.T):
        raise ValueError("invalid covariance")
    W=np.linalg.solve(C,X)
    cov=np.linalg.inv(X.T@W)
    a,b,_=cov@(X.T@np.linalg.solve(C,y))
    if a>=0:
        raise ValueError("not a concave peak")
    vertex=float(np.mean(x)-b/(2*a))
    if not x[0]<=vertex<=x[-1]:
        raise ValueError("peak outside fitted window")
    gradient=np.array([b/(2*a*a),-1/(2*a),0.])
    return vertex,float(np.sqrt(gradient@cov@gradient))

def peak_window_sensitivity(ell,power,sigma,center,half_widths):
    """Vary local fitting window to expose broadband sensitivity."""
    results=[]
    for width in half_widths:
        if width<=0:
            raise ValueError("positive window required")
        mask=np.abs(np.asarray(ell)-center)<=width
        results.append((float(width),*quadratic_peak(np.asarray(ell)[mask],
            np.asarray(power)[mask],np.asarray(sigma)[mask])))
    return results
