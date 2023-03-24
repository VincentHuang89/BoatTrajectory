# 等相位面作图

#%%
import numpy as np
from matplotlib import pyplot as plt
from math import pi,radians,tanh,atan,degrees,asin
from numpy import cos,sin
import pandas as pd
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
def F(k,N):     #公式(36)
    return (k*tanh(k))**0.5/N

def F1(k,N):    #对公式(36)求导
    return (k*tanh(k))**0.5*(0.5*k*(1 - tanh(k)**2) + 0.5*tanh(k))/(k*tanh(k))/N
#eta**2 = k**2-zeta**2

def y(a,k,N):   #公式(19)
    return -a*F1(k,N)/k/(F(k,N)-k*F1(k,N))*(k**2-F(k,N)**2)**0.5

def x(a,k,N):   #公式(20)
    return a*(k-F(k,N)*F1(k,N))/(k*(F(k,N)-k*F1(k,N)))

N=1.001

plt.figure(dpi=300,figsize=(12, 8))
for a in range(1,5):
    X=[]
    Y=[]
    Y1=[]
    f=[]
    f1=[]
    for k in np.arange(0.4,30,0.1):
        
        f.append(F(k,N))
        f1.append(F1(k,N))
        X.append(x(a,k,N))
        Y.append((y(a,k,N)))
        Y1.append(-(y(a,k,N)))
 

    plt.plot(X, (Y))
    plt.plot(X, (Y1))
plt.title('Froude Number: '+str(N))
#plt.xlim(-1,50)
#plt.ylim(-50,50)

#存在的问题，当N小于1时候，文章中fig2的图对不上
#计算散波夹角alpha
k=(Y1[0]-Y1[1])/(X[0]-X[1])
alpha=degrees(atan(k))

#%计算包络线的夹角r
k1=(0-Y1[5])/(0-X[5])
r=degrees(atan(k1))
plt.plot(np.arange(0,50,1),list(-k1*np.arange(0,50,1)),color='black')
plt.show()
# %%分析随着N变化，alpha和r 角度的变化


Alpha=[]
R=[]
N_delta=0.001
for N in np.arange(1+N_delta,3,N_delta):
    a=1
    X=[]
    Y=[]
    Y1=[]
    f=[]
    f1=[]
    for k in np.arange(1,30,0.1):
        f.append(F(k,N))
        f1.append(F1(k,N))
        X.append(x(a,k,N))
        Y.append((y(a,k,N)))
        Y1.append(-(y(a,k,N)))

    k=(Y1[0]-Y1[1])/(X[0]-X[1])
    Alpha.append(degrees(atan(k)))
#计算包络线的夹角r
    k1=(0-Y1[5])/(0-X[5])   #wavenumber=1.8 
    #k1=(0-Y1[5])/(0-X[5])   #wavenumber=1.5 
    R.append(degrees(atan(k1)))

Res1=pd.DataFrame()
Res1['FroudeNumber']=list(np.arange(1+N_delta,3,N_delta))
Res1['Alpha']=Alpha
Res1['R']=R
Res1['Perc']=Res1['R']/Res1['Alpha']
Res1.to_excel('Alpha_R.xlsx')

X=list(Res1['FroudeNumber'])
Y=list(Res1['Alpha'])
Y1=list(Res1['R'])
plt.figure(dpi=600)
plt.xlabel('Froude Number')
plt.ylabel('Degree')
plt.plot(X,Y,label=r'$\alpha$')
plt.plot(X,Y1,label=r'$\gamma$')
plt.legend()

#%%
fig=plt.figure(dpi=300,figsize=(8,4.5))
ax=fig.add_subplot(111)
ax.set_xlabel(r'$\alpha$',fontsize=13)
ax.set_ylabel(r'$\gamma$',fontsize=13,rotation=0)
ax.plot(Y,Y1)
plt.savefig('/home/huangwj/DAS/BoatTrajectory/Paperfig/a_r.pdf',bbox_inches='tight')

#%% 分析随着N变化，alpha和r 角度的变化（r角度重新分析）
