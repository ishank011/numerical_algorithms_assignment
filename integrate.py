import re
import math
import matplotlib.pyplot as plt
import numpy as np

fn="x**2-x+1"

class Integrate():
    global fn
    
    def solve(self,fn,lower,upper,n,method):
        f,a,b=self.func(fn),float(lower),float(upper)
        sum=f(a)+f(b)
        xarray=[a+(b-a)*i/n for i in range(n+1)]
        plotxarray=list(np.arange(a-20,b+20,0.1))
        plotyarray=[f(i) for i in plotxarray]
        maxy=max(plotyarray)
        miny=min(plotyarray)
        farray=[f(a)]
        
        if method=="trapezoid":
            for i in range(1,n):
                y0=a+(b-a)*(i-1)/n
                y=a+(b-a)*i/n
                c=f(y)
                
                farray.append(c)
                self.graph(y0,y,f(y0),c)
                sum+=2*c
            sum=((b-a)/(2*n))*sum
            
            axisx=np.arange(a-20,b+20,1)
            axisy=np.arange(miny/10-5,maxy/10+5)
            plt.plot(axisx,[0]*len(axisx),'k-')
            plt.plot([0]*len(axisy),axisy,'k-')
            
            plt.plot(plotxarray,plotyarray,'r-',linewidth=1)
            plt.axis([a-20,b+20,miny/10-5,maxy/10+5])
            plt.plot([a,b],[0,0],'ko',markersize=8)
            plt.annotate('(%d,0)'%(a),xy=(a,-0.1), xytext=(a,-0.1))
            plt.annotate('(%d,0)'%(b),xy=(b,-0.1), xytext=(b,-0.1))
            plt.annotate('0',xy=(-0.2,-0.2),xytext=(-0.2,-0.2))
            plt.show()
        elif method=="simpson":
            for i in range(1,n):
                y=a+(b-a)*i/n
                c=f(y)
                farray.append(c)
                if i%2==0:
                    sum+=2*c
                else:
                    sum+=4*c
            sum=((b-a)/(3*n))*sum
        
        farray.append(f(b))
        
        return sum
            
            
    
    
        
    def func(self,fn):
        
        reglog=r"(?P<start>.*)\bl(og|n)\(x\)\B(?P<end>.*)"
        reslog=re.search(reglog,fn)
        regsin=r"(?P<start>.*)\bsin\(x\)\B(?P<end>.*)"
        ressin=re.search(regsin,fn)
        regcos=r"(?P<start>.*)\bcos\(x\)\B(?P<end>.*)"
        rescos=re.search(regcos,fn)
        
        
        if reslog!=None:
            if reslog.group('start')=='' and reslog.group('end')=='':
                f=math.log
                print('HURRAY')
            else:
                raise 'Function Type Error'
                
        elif ressin!=None:
            if ressin.group('start')=='' and ressin.group('end')=='':
                f=math.sin
            else:
                raise 'Function Type Error'
        elif rescos!=None:
            if rescos.group('start')=='' and rescos.group('end')=='':
                f=math.cos
            else:
                raise 'Function Type Error'
        else:
            
            f=self.pol
        return f
        
    def polynomial(fn):
        reg=r"([+-]?(\d)*(\*)?x?(\*\*|\^)?(\d)*)"
        a=re.findall(reg,fn)
        pol=[[],[]]
        for i in range(len(a)):
            minus=0
            b=a[i][0]
            
            if b!='':
                if b[0]=='+':
                    b=b[1:]
                elif b[0]=='-':
                    minus=1
                    b=b[1:]
                
                if 'x' not in b:
                    coeff=float(b)
                    power=0
                else:
                    c=re.findall(r"(?P<coeff>\d+)?\*?x(?P<order>(\*\*|\^)(\d+))?",b)
                    
                    if c[0][0]=='':
                        coeff=1.0
                    else:
                        coeff=float(c[0][0])
                    
                    if c[0][-1]=='':
                        power=1.0
                        
                    else:
                        power=float(c[0][-1])
                    
                if minus==1:
                    coeff=-coeff
                
                pol[0].append(power)
                pol[1].append(coeff)
    
        return pol
           
        
        
                
           
    def pol(self,x,power=polynomial(fn)[0],coeff=polynomial(fn)[1]):
        ans=sum([coeff[i]*x**(power[i]) for i in range(len(power))])
        return ans
                
    def graph(self,x1,x2,y1,y2):
        plt.plot([x1,x1,x2,x2],[0,y1,y2,0],'b-',linewidth=2)
        plt.fill_between([x1,x1,x2,x2],[0,y1,y2,0],color='blue',alpha='0.1')

f=fn
n=1000
lower=-4
upper=8
igr=Integrate()
for method in ["trapezoid","simpson"]:
    solution=igr.solve(f,lower,upper,n,method)
    print(solution)
