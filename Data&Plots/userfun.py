from scipy.stats import uniform
from scipy.stats import vonmises
import numpy as np
def sd2k (S):
    R = np.exp(-S**2/2)
    K = 1./(R**3 - 4 * R**2 + 3 * R)
#     K[R < 0.85] = -0.4 + 1.39 * R[R < 0.85] + 0.43/(1 - R[R < 0.85]);
#     K[R < 0.53] = 2 * R[R < 0.53] + R[R < 0.53]**3 + (5 * R[R < 0.53]**5)/6;
    if R<0.53:
        K = 2*R+R**3+(5*R**5)/6
    elif R<0.85:
        K = -0.4+1.39*R+0.43/(1-R)
    return K
def deg2k(sd):
    k = sd2k(np.deg2rad(sd))
    return k
def smi(data,g,B1,B2,sd,mu):
    data = np.deg2rad(data)
    rv = vonmises(deg2k(sd),loc = np.deg2rad(mu))
    rv1 = vonmises(deg2k(sd),loc = np.pi/2)
    rv2 = vonmises(deg2k(sd),loc = -np.pi/2)
    gv = uniform(loc = data[0],scale = np.rad2deg(data[-1]-data[0]))
    y = (1-g-B1-B2)*rv.pdf(data)/360*2*np.pi + (g)*gv.pdf(data) + B1*rv1.pdf(data)/360*2*np.pi + B2*rv2.pdf(data)/360*2*np.pi
    return y