from sys import path; path.append('/home/alberto/esAsi/kolmogorovTools')
from datetime import datetime
import kolmogorovTools as kg
import matplotlib.pyplot as plt
from numpy import linspace, meshgrid, sqrt, ones, savetxt, pi, real, imag
from scipy.sparse.linalg import eigs


def fieldX(x, y, alpha, a):

    fieldX = alpha*x*(1 - (x**2 + y**2)) - y*(1 + alpha*a*(x**2 + y**2))
    return fieldX


def fieldY(x, y, alpha, a):

    fieldY = alpha*y*(1 - (x**2 + y**2)) + x*(1 + alpha*a*(x**2 + y**2))
    return fieldY


def plotSpectrum(eigenVals):
    
    plt.axhline(y=0.0, color='k', linestyle='--', alpha = 0.65, zorder=1)
    plt.axvline(x=0.0, color='k', linestyle='--', alpha = 0.65, zorder=1)
    
    plt.scatter(real(eigenVals), imag(eigenVals), c='C9', marker="o", s=150, zorder=2)        
    plt.xlabel(r'$Re(\lambda)$', size = 26); plt.ylabel(r'$Im(\lambda)$', size = 26)
    plt.tick_params(axis='both', labelsize=22)

    plt.show()


def main():
    
    startTime = datetime.now()
    
    # We define the noise and the parameters
    D = 0.04
    alpha, a, gx1, gx2, gy1, gy2  = [1, -0.3, sqrt(2*D), 0, 0, sqrt(2*D)];
    
    # We define the coordinates of our domain (yp = yPlus, ym= yMinus) - same for xp and xm 
    yp = 2; ym = -2
    xp = 2; xm = -2
    y = linspace(ym, yp, 251); M = len(y); dy = y[1] - y[0]  # Choose how many points you want for your discretisation
    x = linspace(xm, xp, 251); N = len(x); dx = x[1] - x[0]

    # We create the grid
    thetaa, sigmaa = meshgrid(x, y, sparse=True)
    # We evaluate the fields x and y
    xField = fieldX(thetaa, sigmaa, alpha, a).reshape(N*M)
    yField = fieldY(thetaa, sigmaa, alpha, a).reshape(N*M)

    '''
    Next, we create the terms composing the Ldagger operator
    Please note: the boolean terms in the functions mean True = implement reflecting Boundary Conditions. 
    Also note that in the gradXX (gradYY) term we insert a vector of ones because the function g(x) (g(y)) is a constant function (additive noise)
    in the Langevin equation we considered in this example.
    Multiplicative noise is possible to consider by substituting the vector of ones by the corresponding function g(x)^2 (g(y)^2) 
    '''

    # nablaX stands for the term fx*partial_x 
    nablaX = kg.gradX(xField/dx, N, M, True);  
    # nablaXX stands for the term gxx*partial_{xx} 
    nablaXX = kg.gradXX(ones(N*M)/(dx**2), N, M, True) 
    # nablaY stands for the term fy*partial_y
    nablaY = kg.gradY(yField/dy, N, M, True); 
    # nablaYY stands for the term gyy*partial_{yy}
    nablaYY = kg.gradYY(ones(N*M)/(dy**2), N, M, True) 
    
    # We build the lDagger operator
    lDagger = nablaX + nablaY + 0.5*((gx1**2)*nablaXX + (gy2**2)*nablaYY)
    # We diagonalise it
    eigenVals, eigenVects = eigs(lDagger, k=23, sigma=-0.16)
    
    plotSpectrum(eigenVals) # Comment this block if you prefer not to plot the computed spectrum

    # Get eigenfunctions of the backwards operator
    strToSave = './eigenVals_%s' % D; savetxt(strToSave, eigenVals)    
    strToSave = './eigenVects_%s' % D; savetxt(strToSave, eigenVects)

    # Get eigenfunctions of the forwards operator
    eigenVals, eigenVects = eigs(lDagger.H, k=5, sigma=0)
    strToSave = './eigenValsFP_%s' % D; savetxt(strToSave, eigenVals)
    strToSave = './eigenVectsFP_%s' % D; savetxt(strToSave, eigenVects)

    print('\tiniTime: %s\n\tendTime: %s' % (startTime, datetime.now()))
    
    
if __name__ == '__main__':
    main()
    
    
