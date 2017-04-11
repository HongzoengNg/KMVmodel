#calcul totalAsset and thetaAsset
from scipy.optimize import fsolve
import math

def st(x):
    x=float(x)
    return math.exp(-(x**2)/2)

def N(x):
    from scipy import integrate
    nx,err=integrate.quad(st,-100,x)
    nx=nx/math.sqrt(2*math.pi)
    return nx

def f(Ve,ThetaE,D,r,t=1):
    Ve=float(Ve)
    ThetaE=float(ThetaE)
    D=float(D)
    r=float(r)
    t=float(t)
    def ff(x):
        x0=float(x[0])
        x1=float(x[1])
        d1=math.log(abs(x0))-math.log(D)+(r+0.5*(x1**2))*t
        d2=math.log(abs(x0))-math.log(D)+(r+0.5*(x1**2))*t-x1*math.sqrt(t)
        return [
            x0*N(d1)-D*math.exp(-r*t)*N(d2)-Ve,
            x0*x1*(N(d1)/Ve)-ThetaE
        ]
    result=fsolve(ff,[Ve,ThetaE])
    return result

