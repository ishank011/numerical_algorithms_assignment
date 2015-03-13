def f(order, coefficients, x):
    s = sum([coefficients[i]*(x**(order-i)) for i in range(order+1)])
    return s

def diff(order, coefficients, x):
    s = sum([coefficients[i]*(order-i)*(x**(order-i-1)) for i in range(order)])
    return s


class Integrate(object):
    def Trapezoid(self, order, coefficients, low, high, interval):
        N = (high - low)/interval
        I = sum([2*f(order, coefficients, low + (i*interval)) for i in range(1, int(N))])
        I += (f(order, coefficients, low) + f(order, coefficients, high))
        return (I*(high - low))/(2*N)

    def Simpson(self, order, coefficients, low, high, interval):
        N = (high - low)/(2*interval)
        I = sum([4*f(order, coefficients, low + (i*interval)) if i%2==1 else 2*f(order, coefficients, low + (i*interval)) for i in range(1, int(2*N))])
        I += (f(order, coefficients, low) + f(order, coefficients, high))
        return (I*(high - low))/(6*N)
    
    def solve(self, order, coefficients, low, high, method, interval = 1e-3):
        if method == 'trapezoid':
            return self.Trapezoid(order, coefficients, low, high, interval)
        elif method == 'simpson':
            return self.Simpson(order, coefficients, low, high, interval)


class LinearSystemSolver(object):
    def Gauss(self, A, b):
        n, x = len(A), []
        for i in range(n):
            A[i].append(b[i])
        for i in range(1, n):
            y = n-1
            while A[i-1][i-1] == 0:
                A[i-1], A[y] = A[y], A[i-1]
                y-=1
            for j in range(i, n):
                z = A[j][i-1]/A[i-1][i-1]
                for k in range(i-1, n+1):
                    A[j][k]-=(z*A[i-1][k])
        for i in range(n-1, -1, -1):
            if A[i][i] == 0:
                return None
            s = sum([A[i][j]*x[n-j-1] for j in range(i+1, n)])
            x.append((A[i][n]-s)/A[i][i])
        return list(reversed(x))

    def GaussJordan(self, A, b):
        n, x = len(A), []
        for i in range(n):
            A[i].append(b[i])
        for i in range(1, n+1):
            y = n-1
            while A[i-1][i-1] == 0:
                A[i-1], A[y] = A[y], A[i-1]
                y-=1
            for j in range(0, n):
                if j!=(i-1):
                    z = A[j][i-1]/A[i-1][i-1]
                    for k in range(i-1, n+1):
                        A[j][k]-=(z*A[i-1][k])
        for i in range(n-1, -1, -1):
            if A[i][i] == 0:
                return None
            x.append(A[i][n]/A[i][i])
        return list(reversed(x))

    import numpy as np
     
    def GaussSiedel(self, A, b, epsilon):
        A = np.array(A)
        b = np.array(b)
        L = np.tril(A)
        U = A - L
        x = np.ones(len(b))
        i = 0
        while i < 500 and not np.allclose(np.dot(A, x), b, epsilon):
            x = np.dot(np.linalg.inv(L), b - np.dot(U, x))
            i+=1
        return x
    
    def solve(self, A, b, method, epsilon = 1e-6):
        if method == 'gauss':
            return self.Gauss(A, b)
        elif method == 'gauss-jordan':
            return self.GaussJordan(A, b)
        elif method == 'gauss-siedel':
            return self.GaussSiedel(A, b, epsilon)


class PolynomialSolver(object):
    def BisectionSearch(self, order, coefficients, low, high, epsilon):
        num = 0
        while (high-low)>=epsilon and num<500:
            mid = (low+high)/2
            if f(coefficients, low)*f(coefficients, mid) > 0:
                low = mid
            else:
                high = mid
            num += 1
        return mid

    def Secant(self, order, coefficients, low, high, epsilon):
        num, x0, x1, x2 = 0, 1, 2, 3
        while abs(f(order, coefficients, x1))>=epsilon and num<500:
            x2, x1, x0 = (x1 - ((x1-x0)*f(order, coefficients, x1))/(f(order, coefficients, x1)-f(order, coefficients, x0))), x2, x1
            num += 1
        return x1

    def secantRF(self, order, coefficients, epsilon):
    num, x0, x1, x2 = 0, 1, 2, 8
    while abs(f(order, coefficients, x1))>=epsilon and num<500:
        x2, x1, x0 = (x1 - ((x1-x0)*f(order, coefficients, x1))/(f(order, coefficients, x1)-f(order, coefficients, x0))), x2, x1
        while f(order, coefficients, x2)*f(order, coefficients, x1) > 0:
            x2, x1, x0 = (x1 - ((x1-x0)*f(order, coefficients, x1))/(f(order, coefficients, x1)-f(order, coefficients, x0))), x2, x1
        num += 1
    return x1
    
    def NewtonRaphson(self, order, coefficients, low, high, epsilon):
        num, x0, x1 = 0, low, high
        while abs(f(order, coefficients, x0))>=epsilon and num<500:
            x1, x0 = (x0 - (f(order, coefficients, x0)/diff(order, coefficients, x0))), x1
            num += 1
        return x0

    def solve(self, order, coefficients, method, low, high, epsilon = 1e-6):
        if method == 'bisection':
            return self.BisectionSearch(order, coefficients, low, high, epsilon)
        elif method == 'secant':
            return self.Secant(order, coefficients, low, high, epsilon)
        elif method == 'secantrf':
            return self.SecantRF(order, coefficients, low, high, epsilon)
        elif method == 'newtonraphson':
            return self.NewtonRaphson(order, coefficients, low, high, epsilon)


class Interpolate:
    
    def solve(self,L,M,method):
        if(method=="newton"):
            return (self.Newton(L,M))
        else:
            return (self.Lagrange(L,M))

    def Lagrange(self,L,M):                                                
        from numpy import array
        from numpy.polynomial import polynomial as P
        n=len(L)                                                           
        w=(-1*L[0],1)                                                      
        for i in range(1,n):
            w=P.polymul(w,(-1*L[i],1))                                    
        result=array([0.0 for i in range(len(w)-1)])                    
        derivative=P.polyder(w)                                             
        for i in range(n):
            result+=(P.polydiv(w,(-1*L[i],1))[0]*M[i])/P.polyval(L[i],derivative)   
        return(list(result))                                                
    def Newton(self,L,M):                                                   
      
        from numpy import array
        from numpy.polynomial import polynomial as P
        n=len(L)                                                            
        mat=[[0.0 for i in range(n)] for j in range(n)]                    
        for i in range(n):                                                 
            mat[i][0]=M[i]
        for i in range(1,n):                                               
            for j in range(n-i):
                mat[j][i]=(mat[j+1][i-1]-mat[j][i-1])/(L[j+i]-L[j])
        result=array((mat[0][0],))                                          
        for i in range(1,n):
            prod=(-1*L[0],1)                                               
                                                                            
            for j in range(1,i):
                prod=P.polymul(prod,(-1*L[j],1))                              
            result=P.polyadd(result,array(prod)*mat[0][i])                  
        return (list(result))                                               


import numpy as np

class LPsolver():

    def solve(self,a,b,c):
        l_b=len(b)
        l_a=len(a)
        arr=np.array(c)
        arr=np.concatenate((arr,np.identity(len(b))),1)
        b=np.array(b)
        x=[0]*l_a

        arr=np.insert(arr,len(a)+l_b,b,1)
        B=np.array([0]*l_b)
        C=np.array(a+[0]*l_b)
        bx=list(range(l_a,l_a+l_b))
        while(1):
            maxpos=[]
            minpos=[]
            for i in range(l_a+l_b):
                maxpos.append(C[i]-np.dot(B,arr[:,i]))

            c=maxpos.index(max(maxpos))

            if maxpos[c]<=0:
                break
            for i in range(l_b):
                if arr[i,c]!=0:
                    minpos.append(arr[i,-1]/arr[i,c])
                else:
                    minpos.append(10000000000)
            q=filter(lambda x:x>0,minpos)
            if q==[]:
                return "unbounded function"
                break
            r=minpos.index(min(q))

            B[r]=C[c]
            bx[r]=c
            arr[r,:]=arr[r,:]/arr[r,c]
            for i in range(l_b):
                if(i!=r):
                    arr[i,:]=arr[i,:]-(arr[i,c]/arr[r,c])*arr[r,:]

        for i in range(len(bx)):
            if bx[i]<l_a:
                x[bx[i]]=arr[i,-1]
        x=list(map(int,x))

        x.append(sum([x[i]*a[i] for i in range(l_a)]))
        x[-1]=(x[-1],)
        return x
