from numpy import ones, zeros
from scipy.sparse import kron, diags, eye
from scipy.sparse.linalg import lsqr, spsolve


def gradX(xVector, N, M, boundary = False):
    
    # Las derivadas centrales
    center4 = diags((1/280.)*ones(N-4),-4) - diags((1/280.)*ones(N-4),4) - diags((4/105.)*ones(N-3),-3) + diags((4/105.)*ones(N-3),3) + diags((1/5.)*ones(N-2),-2) - diags((1/5.)*ones(N-2),2) - diags((4/5.)*ones(N-1),-1) + diags((4/5.)*ones(N-1),1)
    vectorUseful = zeros(N); vectorUseful[4:N-4] = 1
    center4 = center4.multiply(vectorUseful.reshape(N,1))

    center3 = -diags((1/60.)*ones(N-3),-3) + diags((1/60.)*ones(N-3),3) + diags((3/20.)*ones(N-2),-2) - diags((3/20.)*ones(N-2),2) - diags((3/4.)*ones(N-1),-1) + diags((3/4.)*ones(N-1),1)
    vectorUseful = zeros(N); vectorUseful[3] = 1; vectorUseful[-4] = 1
    center3 = center3.multiply(vectorUseful.reshape(N,1))
    
    center2 = diags((1/12.)*ones(N-2),-2) - diags((1/12.)*ones(N-2),2) - diags((2/3.)*ones(N-1),-1) + diags((2/3.)*ones(N-1),1)
    vectorUseful = zeros(N); vectorUseful[2] = 1; vectorUseful[-3] = 1 
    center2 = center2.multiply(vectorUseful.reshape(N,1))
    
    center1 = -diags((1/2)*ones(N-1),-1) + diags((1/2)*ones(N-1),1)
    vectorUseful = zeros(N); vectorUseful[1] = 1; vectorUseful[-2] = 1 
    center1 = center1.multiply(vectorUseful.reshape(N,1))
    
    # Bloque de arriba y de abajo (forward/backward)
    
    upTerm = -diags((49/20.)*ones(N),0) + diags(6*ones(N-1),1) - diags((15./2)*ones(N-2),2) + diags((20/3.)*ones(N-3),3) - diags((15/4.)*ones(N-4),4) + diags((6/5.)*ones(N-5),5) - diags((1/6.)*ones(N-6),6)
    vectorUseful = zeros(N); vectorUseful[0] = 1
    upTerm = upTerm.multiply(vectorUseful.reshape(N,1))
    downTerm = diags((49/20.)*ones(N),0) - diags(6*ones(N-1),-1) + diags((15./2)*ones(N-2),-2) - diags((20/3.)*ones(N-3),-3) + diags((15/4.)*ones(N-4),-4) - diags((6/5.)*ones(N-5),-5) + diags((1/6.)*ones(N-6),-6)
    vectorUseful = zeros(N); vectorUseful[-1] = 1
    downTerm = downTerm.multiply(vectorUseful.reshape(N,1))   
   
    if boundary: 
        xTerm = center4 + center3 + center2 + center1
    else:
        xTerm = center4 + center3 + center2 + center1 + upTerm + downTerm
    
    xTerm = kron(eye(M), xTerm)    
    xTerm = xTerm.multiply(xVector.reshape(N*M,1))    
    
    return xTerm


def gradY(yVector, N, M, boundary = False):
    
    # Las derivadas centrales
    center4 = diags((1/280.)*ones(N*M-4*N),-4*N) - diags((1/280.)*ones(N*M-4*N),4*N) -diags((4/105.)*ones(N*M-3*N),-3*N) + diags((4/105.)*ones(N*M-3*N),3*N) + diags((1/5.)*ones(N*M-2*N),-2*N) - diags((1/5.)*ones(N*M-2*N),2*N) - diags((4/5.)*ones(N*M-N),-N) + diags((4/5.)*ones(N*M-N),N)
    vectorUseful = zeros(N*M); vectorUseful[4*N:N*(M-4)] = 1
    center4 = center4.multiply(vectorUseful.reshape(N*M,1))    
    
    center3 = -diags((1/60.)*ones(N*M-3*N),-3*N) + diags((1/60.)*ones(N*M-3*N),3*N) + diags((3/20.)*ones(N*M-2*N),-2*N) - diags((3/20.)*ones(N*M-2*N),2*N) - diags((3/4.)*ones(N*M-N),-N) + diags((3/4.)*ones(N*M-N),N)
    vectorUseful = zeros(N*M); vectorUseful[3*N:4*N] = 1; vectorUseful[N*(M-4):N*(M-3)] = 1 
    center3 = center3.multiply(vectorUseful.reshape(N*M,1))
    
    center2 = diags((1/12.)*ones(N*M-2*N),-2*N) - diags((1/12.)*ones(N*M-2*N),2*N) - diags((2/3.)*ones(N*M-N),-N) + diags((2/3.)*ones(N*M-N),N)
    vectorUseful = zeros(N*M); vectorUseful[2*N:3*N] = 1; vectorUseful[N*(M-3):N*(M-2)] = 1 
    center2 = center2.multiply(vectorUseful.reshape(N*M,1))
    
    center1 = -diags((1/2.)*ones(N*M-N),-N) + diags((1/2.)*ones(N*M-N),N)
    vectorUseful = zeros(N*M); vectorUseful[N:2*N] = 1; vectorUseful[N*(M-2):N*(M-1)] = 1 
    center1 = center1.multiply(vectorUseful.reshape(N*M,1))
    
    # Bloque de arriba y de abajo (forward/backward)
    upTerm = -diags((49/20.)*ones(N*M),0) + diags(6*ones(N*M-N),N) - diags((15/2)*ones(N*M-2*N),2*N) + diags((20/3.)*ones(N*M-3*N),3*N) - diags((15/4.)*ones(N*M-4*N),4*N) + diags((6/5.)*ones(N*M-5*N),5*N) - diags((1/6.)*ones(N*M-6*N),6*N)
    vectorUseful = zeros(N*M); vectorUseful[0:N] = 1
    upTerm = upTerm.multiply(vectorUseful.reshape(N*M,1))
    downTerm = diags((49/20.)*ones(N*M),0) - diags(6*ones(N*M-N),-N) + diags((15/2)*ones(N*M-2*N),-2*N) - diags((20/3.)*ones(N*M-3*N),-3*N) + diags((15/4.)*ones(N*M-4*N),-4*N) - diags((6/5.)*ones(N*M-5*N),-5*N) + diags((1/6.)*ones(N*M-6*N),-6*N)
    vectorUseful = zeros(N*M); vectorUseful[N*M-N:] = 1
    downTerm = downTerm.multiply(vectorUseful.reshape(N*M,1))   
    
    if boundary: 
        yTerm = center4 + center3 + center2 + center1
    else:
        yTerm = center4 + center3 + center2 + center1 + upTerm + downTerm
        
    yTerm = yTerm.multiply(yVector.reshape(N*M,1))

    return yTerm


def gradXX(xVector, N, M, boundary = False):
    
    # Las derivadas centrales
    center4 = -diags((1/560.)*ones(N-4),-4) - diags((1/560.)*ones(N-4),4) + diags((8/315.)*ones(N-3),-3) + diags((8/315.)*ones(N-3),3) - diags((1/5.)*ones(N-2),-2) - diags((1/5.)*ones(N-2),2) + diags((8/5.)*ones(N-1),-1) + diags((8/5.)*ones(N-1),1) - diags((205/72.)*ones(N),0)
    vectorUseful = zeros(N); vectorUseful[4:N-4] = 1
    center4 = center4.multiply(vectorUseful.reshape(N,1))
    
    center3 = diags((1/90.)*ones(N-3),-3) + diags((1/90.)*ones(N-3),3) - diags((3/20.)*ones(N-2),-2) - diags((3/20.)*ones(N-2),2) + diags((3/2.)*ones(N-1),-1) + diags((3/2.)*ones(N-1),1) - diags((49/18.)*ones(N),0)
    vectorUseful = zeros(N); vectorUseful[3] = 1; vectorUseful[-4] = 1
    center3 = center3.multiply(vectorUseful.reshape(N,1))
    
    center2 = -diags((1/12.)*ones(N-2),-2) - diags((1/12.)*ones(N-2),2) + diags((4/3.)*ones(N-1),-1) + diags((4/3.)*ones(N-1),1) - diags((5/2.)*ones(N),0)
    vectorUseful = zeros(N); vectorUseful[2] = 1; vectorUseful[-3] = 1 
    center2 = center2.multiply(vectorUseful.reshape(N,1))

    center1 = diags(ones(N-1),-1) + diags(ones(N-1),1) - diags(2*ones(N),0)
    vectorUseful = zeros(N); vectorUseful[1] = 1; vectorUseful[-2] = 1 
    center1 = center1.multiply(vectorUseful.reshape(N,1))
    
    center0 = -diags(2*ones(N),0) + diags(2*ones(N-1),-1) + diags(2*ones(N-1),1)
    vectorUseful = zeros(N); vectorUseful[0] = 1; vectorUseful[-1] = 1;
    center0 = center0.multiply(vectorUseful.reshape(N,1))
    
    # Bloque de arriba y de abajo (forward/backward)
    upTerm = diags((469/90.)*ones(N),0) - diags((223/10.)*ones(N-1),1) + diags((879/20.)*ones(N-2),2) - diags((949/18.)*ones(N-3),3) + diags((41)*ones(N-4),4) - diags((201/10.)*ones(N-5),5) + diags((1019/180)*ones(N-6),6) - diags((7/10)*ones(N-7),7)
    vectorUseful = zeros(N); vectorUseful[0] = 1
    upTerm = upTerm.multiply(vectorUseful.reshape(N,1))
    downTerm = diags((469/90.)*ones(N),0) - diags((223/10.)*ones(N-1),-1) + diags((879/20.)*ones(N-2),-2) - diags((949/18.)*ones(N-3),-3) + diags((41)*ones(N-4),-4) - diags((201/10.)*ones(N-5),-5) + diags((1019/180)*ones(N-6),-6) - diags((7/10)*ones(N-7),-7)
    vectorUseful = zeros(N); vectorUseful[-1] = 1
    downTerm = downTerm.multiply(vectorUseful.reshape(N,1))   

    if boundary: 
        xTerm = center4 + center3 + center2 + center1 + center0
    else:
        xTerm = center4 + center3 + center2 + center1 + upTerm + downTerm

    xTerm = kron(eye(M), xTerm)    
    xTerm = xTerm.multiply(xVector.reshape(N*M,1))    

    return xTerm


def gradYY(yVector, N, M, boundary = False):
    
    center4 = -diags((1/560.)*ones(N*M-4*N),-4*N) - diags((1/560.)*ones(N*M-4*N),4*N) + diags((8/315.)*ones(N*M-3*N),-3*N) + diags((8/315.)*ones(N*M-3*N),3*N) - diags((1/5.)*ones(N*M-2*N),-2*N) - diags((1/5.)*ones(N*M-2*N),2*N) + diags((8/5.)*ones(N*M-N),-N) + diags((8/5.)*ones(N*M-N),N) - diags((205/72.)*ones(N*M),0)
    vectorUseful = zeros(N*M); vectorUseful[4*N:N*(M-4)] = 1
    center4 = center4.multiply(vectorUseful.reshape(N*M,1))

    center3 = diags((1/90.)*ones(N*M-3*N),-3*N) + diags((1/90.)*ones(N*M-3*N),3*N) - diags((3/20.)*ones(N*M-2*N),-2*N) - diags((3/20.)*ones(N*M-2*N),2*N) + diags((3/2.)*ones(N*M-N),-N) + diags((3/2.)*ones(N*M-N),N) - diags((49/18.)*ones(N*M),0)
    vectorUseful = zeros(N*M); vectorUseful[3*N:4*N] = 1; vectorUseful[N*(M-4):N*(M-3)] = 1
    center3 = center3.multiply(vectorUseful.reshape(N*M,1))

    center2 = -diags((1/12.)*ones(N*M-2*N),-2*N) - diags((1/12.)*ones(N*M-2*N),2*N) + diags((4/3.)*ones(N*M-N),-N) + diags((4/3.)*ones(N*M-N),N) - diags((5/2.)*ones(N*M),0)
    vectorUseful = zeros(N*M); vectorUseful[2*N:3*N] = 1; vectorUseful[N*(M-3):N*(M-2)] = 1 
    center2 = center2.multiply(vectorUseful.reshape(N*M,1))
    
    center1 = diags(ones(N*M-N),-N) + diags(ones(N*M-N),N) - diags(2*ones(N*M),0)
    vectorUseful = zeros(N*M); vectorUseful[N:2*N] = 1; vectorUseful[N*(M-2):N*(M-1)] = 1 
    center1 = center1.multiply(vectorUseful.reshape(N*M,1))    
    
    center0 = -diags(2*ones(N*M),0) + diags(2*ones(N*M-N),-N) + diags(2*ones(N*M-N),N)
    vectorUseful = zeros(N*M); vectorUseful[0:N] = 1; vectorUseful[N*(M-1):] = 1;
    center0 = center0.multiply(vectorUseful.reshape(N*M,1))
    
    
    upTerm = diags((469/90)*ones(N*M),0) - diags((223/10)*ones(N*M-N),N) + diags((879/20)*ones(N*M-2*N),2*N) - diags((949/18)*ones(N*M-3*N),3*N) + diags((41)*ones(N*M-4*N),4*N) - diags((201/10)*ones(N*M-5*N),5*N) + diags((1019/180)*ones(N*M-6*N),6*N) - diags((7/10)*ones(N*M-7*N),7*N)
    vectorUseful = zeros(N*M); vectorUseful[0:N] = 1
    upTerm = upTerm.multiply(vectorUseful.reshape(N*M,1))
    downTerm = diags((469/90)*ones(N*M),0) - diags((223/10)*ones(N*M-N),-N) + diags((879/20)*ones(N*M-2*N),-2*N) - diags((949/18)*ones(N*M-3*N),-3*N) + diags((41)*ones(N*M-4*N),-4*N) - diags((201/10)*ones(N*M-5*N),-5*N) + diags((1019/180)*ones(N*M-6*N),-6*N) - diags((7/10)*ones(N*M-7*N),-7*N)
    vectorUseful = zeros(N*M); vectorUseful[N*M-N:] = 1
    downTerm = downTerm.multiply(vectorUseful.reshape(N*M,1))   

    if boundary: 
        yTerm = center4 + center3 + center2 + center1 + center0
    else:
        yTerm = center4 + center3 + center2 + center1 + upTerm + downTerm
        
    yTerm = yTerm.multiply(yVector.reshape(N*M,1))

    return yTerm


    
