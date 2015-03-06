def f(order, coefficients, x):
    s = sum([coefficients[i]*(x**(order-i)) for i in range(order+1)])
    return s

def diff(order, coefficients, x):
    s = sum([coefficients[i]*(order-i)*(x**(order-i-1)) for i in range(order)])
    return s

class Integrate(object):
    def Trapezoid(self, order, coefficients, low, high, interval = 1e-3):
        N = (high - low)/interval
        I = sum([2*f(order, coefficients, low + (i*interval)) for i in range(1, int(N))])
        I += (f(order, coefficients, low) + f(order, coefficients, high))
        return (I*(high - low))/(2*N)

    def Simpson(self, order, coefficients, low, high, interval = 1e-3):
        N = (high - low)/(2*interval)
        I = sum([4*f(order, coefficients, low + (i*interval)) if i%2==1 else 2*f(order, coefficients, low + (i*interval)) for i in range(1, int(2*N))])
        I += (f(order, coefficients, low) + f(order, coefficients, high))
        return (I*(high - low))/(6*N)
    
    def solve(self, order, coefficients, low, high, method, interval = 1e-3):
        if method == 'trapezoid':
            return self.Trapezoid(order, coefficients, low, high, interval)
        elif method == 'simpson':
            return self.Simpson(order, coefficients, low, high, interval)

#Gauss-Siedel iterations not done yet
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
    
    def solve(self, A, b, method):
        if method == 'gauss':
            return self.Gauss(A, b)
        elif method == 'gauss-jordan':
            return self.GaussJordan(A, b)
        elif method == 'gauss-siedel':
            return self.GaussSiedel(A, b)


#Regula Falsi method not done yet
class PolynomialSolver(object):
    def BisectionSearch(order, coefficients, low, high, epsilon = 1e-6):
        num = 0
        while (high-low)>=epsilon and num<500:
            mid = (low+high)/2
            if f(coefficients, low)*f(coefficients, mid) > 0:
                low = mid
            else:
                high = mid
            num += 1
        return mid

    def Secant(order, coefficients, low, high, epsilon = 1e-6):
        num, x0, x1, x2 = 0, 1, 2, 3
        while abs(f(order, coefficients, x1))>=epsilon and num<500:
            x2, x1, x0 = (x1 - ((x1-x0)*f(order, coefficients, x1))/(f(order, coefficients, x1)-f(order, coefficients, x0))), x2, x1
            num += 1
        return x1
    
    def NewtonRaphson(order, coefficients, low, high, epsilon = 1e-6):
        num, x0, x1 = 0, low, high
        while abs(f(order, coefficients, x0))>=epsilon and num<500:
            x1, x0 = (x0 - (f(order, coefficients, x0)/diff(order, coefficients, x0))), x1
            num += 1
        return x0

    def solve(self, order, coefficients, method, low, high, epsilon = 1e-6):
        if method == 'bisection':
            return BisectionSearch(order, coefficients, low, high, epsilon = 0.0001)
        elif method == 'secant':
            return Secant(order, coefficients, low, high, epsilon = 0.0001)
        elif method == 'secantrf':
            return SecantRF(order, coefficients, low, high, epsilon = 0.0001)
        elif method == 'newtonraphson':
            return NewtonRaphson(order, coefficients, low, high, epsilon = 0.0001)


#Newton's method not done
class Interpolate(object):
    from functools import reduce
    
    def Lagarange(x_values, y_values, x):
        n = len(x_values)
        def basis(j):
            p = [(x - x_values[k])/(x_values[j] - x_values[k]) for k in range(n) if k!=j]
            return reduce(lambda a, b: a*b, p)
        return sum(basis(j)*y_values[j] for j in range(n)

    def solve(self, x_values, y_values):
        if method == 'newton':
            return self.Newton(x_values, y_values)
        elif method == 'lagarange':
            return self.Lagarange(x_values, y_values)
