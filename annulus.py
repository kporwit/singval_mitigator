import numpy as np

def annulus(R1, R2):
    #random throw of theta angle and length of the vector
    deg = np.random.uniform(0.0, np.multiply(2, np.pi))
    r = np.random.uniform(R1, R2)
    #convertion to Cartesian coordinates
    x = np.multiply(r, np.cos(deg))
    y = np.multiply(r, np.sin(deg))
    #return complex number from Cartesian coordinates
    return complex(x,y)

