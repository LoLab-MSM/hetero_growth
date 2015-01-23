import numpy as np
from scipy.special import erf

def pdf(xpoints, mean, sdev, skew):
    delta2 = 0.5*np.pi*np.abs(skew)**(2./3) / ( np.abs(skew)**(2./3) + (0.5*(4-np.pi))**(2./3) )
    delta = np.sign(skew)*np.sqrt(delta2)
    alpha = delta/np.sqrt(1-delta2)
    omega = sdev/np.sqrt(1-2.*delta2/np.pi)
    xi = mean-omega*delta*np.sqrt(2/np.pi)
    return np.array([1./omega/np.sqrt(2.*np.pi)*np.exp(-0.5*((x-xi)/omega)**2)*(1.+erf(alpha*(x-xi)/np.sqrt(2.)/omega)) for x in xpoints])