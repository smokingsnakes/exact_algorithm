# -*- coding: utf-8 -*-
import numpy as np
import scipy
import scipy.integrate as integrate



S=100
X=100
r=0.2
v=0.9
theta = 0.3
rho= -0.2
k=2
sigma =0.9
t=0
tau=1


def phi_heston(a, v0, v_t, d):

    gamma_a = np.sqrt(k**2 - 2 * sigma**2 * 1j*a)
    gammadt = gamma_a * (tau - t)
    sqrtv0vt = np.sqrt(v0*v_t)
    delta = -k * (tau-t)

    part1 = (gamma_a * np.exp(-(gamma_a - k)/2 * (tau-t)) * 
             (1 - np.exp(delta)))/(k * (1- np.exp(- gammadt)))

    part2 = np.exp((v0+v_t)/(sigma**2) * 
                   ( (k * (1 + np.exp(delta)))/(1-np.exp(delta))) -
                   (gamma_a * (1 + np.exp(- gammadt)))/(1 - np.exp(- gammadt)))
    

    part3 = scipy.special.jn( 0.5*d - 1, ((4 * gamma_a * sqrtv0vt)/(sigma**2) *
            np.exp(- gammadt/2)/(1 - np.exp(- gammadt)))) / scipy.special.jn(
                                    0.5*d - 1, 
                                    ((4 * k * sqrtv0vt)/(sigma**2) * 
                             (np.exp(delta/2))/(1 - np.exp(delta))))
    
    return part1 + part2 + part3


def intv(n, cf, v0, v_t, d):
    
    U = np.random.normal(0,1,n)
    
    def integrand(phi = phi_heston, v0 = v0, v_t = v_t, d = d, x = U):
        def f2(u,x):
            func = np.imag(phi(u, v0, v_t, d) * np.exp(-1j * u * x)) / u
            return func
        return f2

    ## integrate to CDF
    
    def F_x(x):
        
        return 0.5 - 1/np.pi * integrate.quad(integrand, 0, np.inf,
                                              args = (x))[0]

    def invcdf(u):
        
        def subcdf(t):
            return F_x(t) - u
        
        return scipy.optimize.newton(subcdf, x0 = 0.1)
    
    
    return invcdf(U)


def heston(S, X, r, v, theta, rho, k, sigma, phi, t = 0, tau = 1):
    
    d1 = (4 * k * theta)/(sigma)**2
    c0 = (sigma**2 * (1 - np.exp(-k*tau)))/(4*k)
    dt = (tau-t)
    ST = None
    
    # sampling V
    lambda1 = (4*k*np.exp(-k*dt)*v)/(sigma**2 * (1-np.exp(-k*dt)))
    vt = c0 * np.random.noncentral_chisquare(size = 1, df = d1, nonc = lambda1)
    
    
    # Sampling int{V}
    int_v = intv(v0 = v, v_t = vt, d = d1, n = 1, cf = phi)
    
    # Sampling int{v}dw
    int_vdw = (1/sigma) * (vt - v - k * theta * dt + k  * int_v)
    
    
    # Sampling S
    if int_v >= 0:
        m = np.log(S) + (r * (tau - t) - (1/2) * int_v + rho * int_vdw)
        std = np.sqrt((1 - rho**2)) * np.sqrt(int_v)
        S = np.exp(m + std * np.random.normal(1))
        v = vt
        ST = S
    else:
        v = vt
        ST = np.array([ST,np.NaN])
    
    Result = ST - X
    Result[Result <= 0] = 0
    call = np.exp(-r * tau) * Result
    
    return(call, Result, ST)





heston(S=S, X=X, r=r, v = v, theta = theta, rho = rho,
       k = k, sigma = sigma,phi = phi_heston)







































