import numpy as np
import cmath

#parameters for morlet functions
sigmas = [1, 3, 6]
thetas = [0, np.pi/6, np.pi/4, np.pi/3, np.pi/2, np.pi*2/3, np.pi*3/4, np.pi*5/6]

# we need to plot two graphs, a real and an imiginary. the morlet function returns a complex number so we need to extract the real and imiginary
# from that function thus the morlet real and imiginary. in order to do that function we need c1 and c2 then we can do that function for each pixel
# the morlet function is a complex function that taat returns a complex number,
# it takes a complex number which is rewritten from e^i etc to cos(something) + i sin(something)

#define the morlet function that return the real part
def morlet_real(x, y, sig, theta, C1, C2):
    # set variables
    cosTheta = np.cos(theta)
    sinTheta = np.sin(theta)
    pie = np.pi
    # one peak morlet wave function with greek letter equal to 4
    exponentOfEInsideBrackets = (pie / (2 * sig)) * ((x * cosTheta) + (y * sinTheta))
    exponentOfEOutsideBrackets = -(x**2 + y**2)/ (2 * sig**2)

    #morlet wave function
    #cmath.rect(r, phi) Return the complex number x with polar coordinates r and phi.
    z = C1 / sig * (cmath.rect(1, exponentOfEInsideBrackets) - C2) * np.exp(exponentOfEOutsideBrackets)
    return z.real


#define the morlet function that return the imaginary part
def morlet_imag(x, y, sig, theta, C1, C2):
    # set variables
    cosTheta = np.cos(theta)
    sinTheta = np.sin(theta)
    pie = np.pi
    # one peak morlet wave function with greek letter equal to 4
    exponentOfEInsideBrackets = (pie / (2 * sig)) * ((x * cosTheta) + (y * sinTheta))
    exponentOfEOutsideBrackets = -(x**2 + y**2)/ (2 * sig**2)

    #morlet wave function
    #cmath.rect(r, phi) Return the complex number x with polar coordinates r and phi.
    z = C1 / sig * (cmath.rect(1, exponentOfEInsideBrackets) - C2) * np.exp(exponentOfEOutsideBrackets)
    return z.imag


#finds the constants c2
def find_c2(xymin, xymax, sig, theta):
    numerator = 0
    denominator = 0
    cosine = np.cos
    cosineTheta = np.cos(theta)
    sineTheta = np.sin(theta)
    pie = np.pi
    for x in range(xymin, xymax+1, 1):
        for y in range( xymin, xymax+1, 1):
            numerator = numerator + (cosine((pie / (2 * sig)) * ((x * cosineTheta) + (y * sineTheta))) * np.exp(-(x**2 + y**2)/(2 * sig**2)))
            denominator = denominator + (np.exp(-(x**2 + y**2)/(2 * sig**2)))

    C2 = numerator/denominator
    return C2


#finds the constant c1
def find_c1(xymin, xymax, sig, theta, C2):
    Z = 0
    pie = np.pi
    cosTheta = np.cos(theta)
    sinTheta = np.sin(theta)
    cosine = np.cos

    for x in range(xymin, xymax+1, 1):
        for y in range( xymin, xymax+1, 1):
            Z = Z + (1 - 2* C2 * cosine(pie/(2*sig) * ((x * cosTheta) + (y * sinTheta))) + C2**2) * np.exp((-(x**2 + y**2)/sig**2))
    C1 = 1/np.sqrt(Z)

    return C1


#plot the morlet function for the real
def morletMatrix_real(xymin, xymax, sig, theta):

    #find c1 and c2
    C2 = find_c2(xymin, xymax, sig, theta)
    C1 = find_c1(xymin, xymax, sig, theta, C2)

    #define grid over which the function should be plotted
    xx, yy = np.meshgrid(np.linspace(xymin, xymax, 33),np.linspace(xymin, xymax, 33))

    # fill a matrix with the morlet function values
    zz= np.zeros(xx.shape)
    for i in range(yy.shape[0]):
        for j in range(xx.shape[0]):
            zz[i,j] = morlet_real(xx[i,j], yy[i,j], sig, theta, C1, C2)

    return zz

# plot morlet function for imiginary
def morletMatrix_imag(xymin, xymax, sig, theta):
    #determine constants
    C2 = find_c2(xymin, xymax, sig, theta)
    C1 = find_c1(xymin, xymax, sig, theta, C2)

    #define grid over which the function should be plotted
    xx = np.meshgrid(np.linspace(xymin, xymax, xymax-xymin+1))
    yy = np.meshgrid(np.linspace(xymin, xymax, xymax-xymin+1))

    # fill a matrix with the morlet function values
    zz= np.zeros(xx.shape)
    for i in range(yy.shape[0]):
        for j in range(xx.shape[0]):
            zz[i,j] = morlet_imag(xx[i,j], yy[i,j], sig, theta, C1, C2)

    return zz

# morletMatrix_imag(33,33,1,0)


