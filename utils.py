import math
import numpy as np

def phase(aa):
    return math.e ** (1j * np.angle(aa))

def mod(aa):
    return np.abs(aa)

def noise_mea(snr, signal):
    noise = np.random.randn(np.shape(signal)[0],np.shape(signal)[1])\
            +1j*np.random.randn(np.shape(signal)[0],np.shape(signal)[1])
    sigma = np.linalg.norm(signal, 'fro')*(10**(-snr/20))/np.linalg.norm(noise, 'fro')
    return signal + sigma*noise


def compute_rdistance(s1,s2):
    return np.linalg.norm(s1 - math.e**(-1j*np.angle(np.trace(np.dot(np.matrix(s1).H,
                                                                   np.matrix(s2)))))*s2)/np.linalg.norm(s1)

