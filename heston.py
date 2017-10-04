# -*- coding: utf-8 -*-
import numpy as np
import scipy
import scipy.integrate as integrate

k=2
sigma=0.2
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
    

    part3 = scipy.special.jn( v = ((4 * gamma_a * sqrtv0vt)/(sigma**2) *
            np.exp(- gammadt/2)/(1 - np.exp(- gammadt))),
                            z = 0.5*d - 1)/ scipy.special.jn( 
                                    v = ((4 * k * sqrtv0vt)/(sigma**2) * 
                             (np.exp(delta/2))/(1 - np.exp(delta))), z = 0.5*d - 1)
    
    return part1 + part2 + part3


def intv(n, cf):
    
    def integrand(x, u, phi):
        
        return np.imag(phi(u) * np.exp(-1j * u * x)) / u

    ## integrate to CDF
    
    def F_x(x):
        
        return 0.5 - 1/np.pi * integrate.quad(integrand(x), 0, np.inf, args = (u))

    def invcdf(u):
        
        def subcdf(t):
            return F_x(t) - u
        
        return scipy.optimize.newton(F_x(t), x0 = 0.1)


