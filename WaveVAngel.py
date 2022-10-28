# 等相位面作图

#%%

from cmath import nan
from turtle import color
import numpy as np
from matplotlib import pyplot as plt
from math import pi,radians,tanh,tan
from numpy import cos,sin
import math
from math import acos,atan,degrees,sqrt,sinh,tanh,asin

from pyparsing import col

def FroudeNum(v,h):
    f=v/sqrt(9.8*h)
    return f


def F(k,N):     #公式(36)
    return (k*tanh(k))**0.5/N

def F1(k,N):    #对公式(36)求导
    return (k*tanh(k))**0.5*(0.5*k*(1 - tanh(k)**2) + 0.5*tanh(k))/(k*tanh(k))/N
#eta**2 = k**2-zeta**2

def y(a,k,N):   #公式(19)
    return -a*F1(k,N)/k/(F(k,N)-k*F1(k,N))*(k**2-F(k,N)**2)**0.5

def x(a,k,N):   #公式(20)
    return a*(k-F(k,N)*F1(k,N))/(k*(F(k,N)-k*F1(k,N)))







#参照论文，根据Froude number与船行波锥形角alpha的曲线，输出特定fr所对应的alpha角
def FRandAngel(fr):
    g=9.8
    H=7
    ALPHA=[]
    FR=[]
    for k in np.arange(0.01,2.5,0.01):
        n=2*k*H/sinh(2*k*H)
        V=(sqrt((tanh(k*H)/(k*H)*(3-n)*g*H)/2))
        Fr=V/(sqrt(g*H))
        FR.append(Fr)
        if V<(sqrt(g*H)):
            alpha= degrees(acos((sqrt(8*(1-n))/(3-n))))
            ALPHA.append(alpha)
    FR.sort()
    ALPHA.sort()

    for Fr in np.arange(1,3,0.01):
            FR.append(Fr)
            alpha=degrees(asin(1/Fr))
            ALPHA.append(alpha)
    FR_interp=np.arange(0,3,0.01)
    ALPHA_interp=np.interp(FR_interp,FR,ALPHA)
    fr_interp=round(fr,2)
    alpha=list(ALPHA_interp)[list(FR_interp).index(fr_interp)]
    return alpha


#输出特定波数a所对应的船行波
def WavePattern(a,n):
    X=[]
    Y=[]
    Y1=[]
    for k in range(1,1000):
        X.append(x(a,k,n))
        Y.append((y(a,k,n)))
        Y1.append(-(y(a,k,n)))
    return X,Y,Y1


def rotate(angle,valuex,valuey):
    rotatex = math.cos(angle)*valuex -math.sin(angle)*valuey
    rotatey = math.cos(angle)*valuey + math.sin(angle)* valuex
    return rotatex,rotatey
def getLen(x1,y1,x2,y2):
    diff_x = (x1-x2)**2
    diff_y = (y1-y2)**2
    length = np.sqrt(diff_x+diff_y)
    return length

    
def ROTATE(angle,X,Y):
    RX=[]
    RY=[]
    for i in range(0,len(X)):
        rotatex,rotatey=rotate(angle,X[i],Y[i])
        RX.append(rotatex)
        RY.append(rotatey)
    return RX,RY

def move(X,Y,angle,dist):
    distX=dist*cos(angle)
    distY=dist*sin(angle)
    X_m=np.array(X)-distX
    Y_m=np.array(Y)-distY
    return X_m,Y_m

def PlotBoatWave(v,h,angle,dist):
    N=FroudeNum(v,h)
    A=5
    plt.figure(figsize=(10,10))
    #plt.xlim(-8,8)
    #plt.ylim(-8,8)
    trajx=np.arange(-10,10,0.1)
    trajy=tan(angle)*trajx
    plt.plot(trajx,trajy,linestyle='--')

    plt.vlines(x=0,ymin=-10,ymax=30,colors='black',linestyle='--')
    plt.legend(["Sailing line","Optical fiber"])
    for a in range(1,A+1):
        X,Y,Y1=WavePattern(a,N)

        RX,RY=ROTATE(angle,X,Y)
        RX1,RY1=ROTATE(angle,X,Y1)

        RMX,RMY=move(RX,RY,angle,dist)
        RMX1,RMY1=move(RX1,RY1,angle,dist)

        lines=plt.plot((RX),(RY))
        clr=lines[0].get_color()
        plt.plot((RX1),(RY1),color=clr)
        lines=plt.plot(RMX,RMY)
        clr=lines[0].get_color()
        plt.plot(RMX1,RMY1,color=clr)


#判断与x=0的交点
def cross(X,Y):
    cp=[]
    for i in range(0,len(X)-1):
        if X[i]*X[i+1]<0:
            cp.append(Y[i]-X[i]*(Y[i]-Y[i+1])/(X[i]-X[i+1]))       
    return cp


def WaveCross(v,h,angle,T): 
    #求解旋转后船行波与x=0（光纤）的多个交点
    A=100    
    N=FroudeNum(v,h)
    DAS=[]
    for t in np.arange(0,T,1):
        dist=v*t/10 #移动的距离 dist=v*t,除以矫正系数
        crossp=[]
        for a in range(1,A):
            X,Y,Y1=WavePattern(a,N)
            RX,RY=ROTATE(angle,X,Y)
            RX1,RY1=ROTATE(angle,X,Y1)
            RMX,RMY=move(RX,RY,angle,dist)
            RMX1,RMY1=move(RX1,RY1,angle,dist)
            cp=cross(RMX,RMY)
            crossp=crossp+cp
            cp=cross(RMX1,RMY1)
            crossp=crossp+cp
        crossp.sort()
        DAS.append(crossp)
    return DAS

def PlotWaveInDAS(crossp):
    for i in range(0,len(crossp)):
        x=[i]*len(crossp[i])
        plt.scatter(x,crossp[i],s=2)
    plt.show



#重构plotWaveInDAS，以a为外循环，dist为内循环
def PlotWaveInDas(v,h,angle,T,A):  
    N=FroudeNum(v,h)
    delta_T=0.01
    plt.figure(dpi=800,figsize=(10,10))
    for a in range(1,A+1):
        crossp1=[]
        crossp2=[]
        for t in np.arange(0,T,delta_T):
            dist=v*t #移动的距离 dist=v*t,除以矫正系数
            X,Y,Y1=WavePattern(a,N)
            RX,RY=ROTATE(angle,X,Y)
            RX1,RY1=ROTATE(angle,X,Y1)
            RMX,RMY=move(RX,RY,angle,dist)
            RMX1,RMY1=move(RX1,RY1,angle,dist)
            cp=cross(RMX,RMY)
            crossp1=crossp1+cp
            cp=cross(RMX1,RMY1)
            crossp2=crossp2+cp
        x1=np.arange(0,len(crossp1),1)
        x1=x1*delta_T
        x2=np.arange(0,len(crossp2),1)
        x2=x2*delta_T

        plt.plot(x1,crossp1,lw=0.5)
        plt.plot(x2,crossp2,lw=0.5)
    plt.savefig('WavePatternInDas.png')