# 等相位面作图

#%%
import numpy as np
from matplotlib import pyplot as plt
from math import pi,radians,tanh
from numpy import cos,sin

theta = radians(60)
theta = np.linspace(-0.5*pi, 0.5*pi, 200)


plt.figure(dpi=300,figsize=(10, 6))
a=1

for a in range(1,8):
    Y=-a*sin(theta)*cos(theta)**2   #公式(32)
    X=a*cos(theta)*(2-cos(theta)**2)  #公式(33)
    plt.plot(X, Y)

plt.show()
#%%

'''
def F(k): #zeta=F(k)
    return k**(1/2)

def F1(k):
    return 0.5*k**(-0.5)    

'''
def F(k):     #公式(36)
    N=1.56
    return (k*tanh(k))**0.5/N

def F1(k):    #对公式(36)求导
    N=1.56
    return (k*tanh(k))**0.5*(0.5*k*(1 - tanh(k)**2) + 0.5*tanh(k))/(k*tanh(k))/N
N=1.56
#eta**2 = k**2-zeta**2

def y(a,k):   #公式(19)
    return -a*F1(k)/k/(F(k)-k*F1(k))*(k**2-F(k)**2)**0.5

def x(a,k):   #公式(20)
    return a*(k-F(k)*F1(k))/(k*(F(k)-k*F1(k)))


plt.figure(dpi=300,figsize=(12, 8))
for a in range(1,6):
    X=[]
    Y=[]
    Y1=[]
    f=[]
    f1=[]
    for k in range(1,10):
        f.append(F(k))
        f1.append(F1(k))
        X.append(x(a,k))
        Y.append((y(a,k)))
        Y1.append(-(y(a,k)))

    plt.plot(X, (Y))
    plt.plot(X, (Y1))
plt.title('Froude Number: '+str(N))
plt.xlim(-1,30)
plt.ylim(-15,15)
plt.show()
#存在的问题，当N小于1时候，文章中fig2的图对不上

# %%
for a in range(4,5):
    X=[]
    Y=[]
    Y1=[]
    f=[]
    f1=[]
    for k in range(1,1000):
        f.append(F(k))
        f1.append(F1(k))
        X.append(x(a,k))
        Y.append((y(a,k)))
        Y1.append(-(y(a,k)))

    plt.plot(X, (Y))
    plt.plot(X, (Y1))
# %%计算斜率
from math import atan,degrees
k=[]
beta=[]
for i in range(1,len(X)):
    k.append((Y[i]-Y[i-1])/(X[i]-X[i-1]))
    beta.append(degrees(atan((Y[i]-Y[i-1])/(X[i]-X[i-1]))))


# %%
