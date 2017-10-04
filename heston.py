# -*- coding: utf-8 -*-
import numpy as np
import scipy


def phi_heston(a, v0, v_t, d):

    gamma_a = np.sqrt(k**2 - 2 8 sigma**2 * 1j*a)
    gammadt = gamma_a * (tau - t)
    sqrtv0vt = np.sqrt(v0*v_t)
    delta = -k * (tau-t)

    part1 = (gamma_a * np.exp(-(gamma_a - k)/2 * (tau-t)) * (1 - exp(delta)))/
        (k * (1- np.exp(- gammadt)))

    part2 = np.exp((v0+v_t)/(sigma**2) * ( (k * (1 + np.exp(delta)))/(1-np.exp(delta))) -
                   (gamma_a * (1 + np.exp(- gammadt)))/(1 - np.exp(- gammadt))))
    

    part3 = scipy.special.jn(((4 * gamma_a * sqrtv0vt)/(sigma**2) *
                                      np.exp(- gammadt/2)/
                                      (1 - np.exp(- gammadt))), 0.5*d - 1) /
            scipy.special.jn(((4 * k * sqrtv0vt)/(sigma**2) * (np.exp(delta/2))/
                                 (1 - np.exp(delta))), 0.5*d - 1)
    
    return part1 + part2 + part3



