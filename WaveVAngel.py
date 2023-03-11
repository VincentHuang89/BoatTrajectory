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
def WavePattern(a,n,k0=1,kmax=25,kdelta=0.01):
    X=[]
    Y=[]
    Y1=[]
    for k in np.arange(k0,kmax,kdelta):
        X.append(x(a,k,n))
        Y.append((y(a,k,n)))
        Y1.append(-(y(a,k,n)))
    #计算波线的角度

    return X,Y,Y1


def WavePattern1(a,n,k0=1,kmax=25,kdelta=0.01):
    X=[]
    Y=[]
    Y1=[]
    for k in np.arange(k0,kmax,kdelta):
        X.append(x(a,k,n))
        Y.append((y(a,k,n)))
        Y1.append(-(y(a,k,n)))
    #计算波线的角度
    alpha=[degrees(atan((Y1[1]-Y1[0])/(X[1]-X[0])))]
    for i in range(1,len(X)):

        alpha.append(degrees(atan((Y1[i]-Y1[i-1])/(X[i]-X[i-1]))))
    return X,Y,Y1,alpha


def rotate(angle,valuex,valuey):
    rotatex = math.cos(radians(angle))*valuex -math.sin(radians(angle))*valuey
    rotatey = math.cos(radians(angle))*valuey + math.sin(radians(angle))* valuex
    rotatex = math.cos(radians(angle))*valuex -math.sin(radians(angle))*valuey
    rotatey = math.cos(radians(angle))*valuey + math.sin(radians(angle))* valuex
    return rotatex,rotatey
def getLen(x1,y1,x2,y2):
    diff_x = (x1-x2)**2
    diff_y = (y1-y2)**2
    length = np.sqrt(diff_x+diff_y)
    return length

    
def ROTATE(angle,X,Y):
    RX=[]
    RY=[]
    angle=angle-90 #重新映射到光纤的角度
    for i in range(0,len(X)):
        rotatex,rotatey=rotate(angle,X[i],Y[i])
        RX.append(rotatex)
        RY.append(rotatey)
    return RX,RY

def move(X,Y,angle,dist):
    angle=angle-90 #重新映射到光纤的角度
    distX=dist*cos(radians(radians(angle)))
    distY=dist*sin(radians(radians(angle)))
    X_m=np.array(X)-distX
    Y_m=np.array(Y)-distY
    return X_m,Y_m

def PlotBoatWave(v,h,angle,dist):
    N=FroudeNum(v,h)
    A=5
    plt.figure(dpi=300,figsize=(10,10))
    #plt.xlim(-8,8)
    #plt.ylim(-8,8)
    trajx=np.arange(-10,10,0.1)
    if angle==90:
        trajy=np.zeros(len(trajx))
    else:
        trajy=tan(radians(radians(angle)))*trajx
    plt.plot(trajx,trajy,linestyle='--')
    plt.xlim((-40,40))
    plt.ylim((-40,40))
    plt.vlines(x=0,ymin=-10,ymax=30,colors='black',linestyle='--')
    plt.xlim([-20,20])
    plt.ylim([-20,20])
    plt.legend(["Sailing line","Optical fiber"])
    for a in range(1,A+1):
        X,Y,Y1=WavePattern(a,N,0.8,20)

        RX,RY=ROTATE(angle,X,Y)
        RX1,RY1=ROTATE(angle,X,Y1)

        RMX,RMY=move(RX,RY,angle,dist)
        RMX1,RMY1=move(RX1,RY1,angle,dist)

        lines=plt.plot((RX),(RY),lw=3)
        clr=lines[0].get_color()
        plt.plot((RX1),(RY1),color=clr,lw=3)
        lines=plt.plot(RMX,RMY,lw=3)
        clr=lines[0].get_color()
        plt.plot(RMX1,RMY1,color=clr,lw=3)


def PlotShipWakePattern(v,h,A,angle,dist):
    '''
        Depict the ship wake

    '''
    N=FroudeNum(v,h)

    plt.figure(dpi=300,figsize=(10,10))

    plt.vlines(x=0,ymin=-10,ymax=30,colors='black',linestyle='--')
    plt.xlim([-20,20])
    plt.ylim([-20,20])
    plt.legend(["Sailing line","Optical fiber"])
    for a in range(1,A+1):
        X,Y,Y1=WavePattern(a,N)
        lines=plt.plot((X),(Y),lw=3)
        clr=lines[0].get_color()
        plt.plot((X),(Y1),color=clr,lw=3)


    plt.figure(dpi=300,figsize=(10,10))
    for a in range(1,A+1):
        X,Y,Y1=WavePattern(a,N)

        RX,RY=ROTATE(angle,X,Y)
        RX1,RY1=ROTATE(angle,X,Y1)

        RMX,RMY=move(RX,RY,angle,dist)
        RMX1,RMY1=move(RX1,RY1,angle,dist)

        lines=plt.plot((RX),(RY),lw=3)
        clr=lines[0].get_color()
        plt.plot((RX1),(RY1),color=clr,lw=3)
        lines=plt.plot(RMX,RMY,lw=3)
        clr=lines[0].get_color()
        plt.plot(RMX1,RMY1,color=clr,lw=3)



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
    delta_T=0.1
    plt.figure(dpi=300,figsize=(10,10))
    for a in range(1,A+1):
        crossp1=[]
        crossp2=[]
        for t in np.arange(0,T,delta_T):
            dist=v*t #移动的距离 dist=v*t,除以矫正系数
            X,Y,Y1=WavePattern(a,N,1,1.5,0.01)
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

        plt.plot(x1,crossp1,lw=3)
        plt.plot(x2,crossp2,lw=3)
    plt.savefig('WavePatternInDas.png')

def PlotSimuWaveInDas(v,h,phy,T,A,WLen_Scale,Wbias,kmin,kmax):  

    UpperBound=1000
    LowerBound=-1000
    N=FroudeNum(v,h)
    delta_T=0.01
    plt.figure(dpi=300,figsize=(6,6))
    for a in range(1,A+1):
        crossp1=[]
        crossp2=[]
        t_start=0
        t_start1=0
        for t in np.arange(0,T,delta_T):
            dist=v*t 
            X,Y,Y1,ALPHA=WavePattern1(a,N,kmin,kmax,0.1)
            alpha=max(ALPHA)
            Attenuation0 = sin(radians(phy - (alpha)))**2
            Attenuation1 = sin(radians(180-phy - (alpha)))**2
            RX,RY=ROTATE(phy,X,Y)
            RX1,RY1=ROTATE(phy,X,Y1)
            RMX,RMY=move(RX,RY,phy,dist)
            RMX1,RMY1=move(RX1,RY1,phy,dist)
            cp=cross(RMX,RMY)
            if len(cp)>0:
                if cp[0]>UpperBound or cp[0]<LowerBound: 
                    cp=[]
            if len(cp) and t_start==0:
                t_start=t
            crossp1=crossp1+cp

            cp=cross(RMX1,RMY1)
            if len(cp)>0:
                if cp[0]>UpperBound or cp[0]<LowerBound: 
                    cp=[]
            if cp and t_start1==0:
                t_start1=t
            crossp2=crossp2+cp
        x1=np.arange(0,len(crossp1),1)
        x1=x1*delta_T+t_start
        x2=np.arange(0,len(crossp2),1)
        x2=x2*delta_T+t_start1

        #Scale the wave length
        x1=WLen_Scale*x1
        x2=WLen_Scale*x2
        crossp1=WLen_Scale*(np.array(crossp1)+Wbias)
        crossp2=WLen_Scale*(np.array(crossp2)+Wbias)
        
        #calculate the sloop speed of divergent wave
        #print('speed1',(crossp1[-1]-crossp1[-2])/(x1[-1]-x1[-2]))
        #print('speed2',(crossp2[-1]-crossp2[-2])/(x2[-1]-x2[-2]))


        
        channel_spacing=1#4.0838
        DownSampleRate=1#50

        crossp1=crossp1/channel_spacing
        crossp2=crossp2/channel_spacing
        x1=x1*DownSampleRate
        x2=x2*DownSampleRate

        plt.plot(x1,crossp1,alpha=Attenuation0,lw=1)
        plt.plot(x2,crossp2,alpha=Attenuation1,lw=1)
        #print(x1,crossp1)
    #plt.ylim([-30,30])
    #plt.xlim([0,6])

    plt.xlabel('Time')
    plt.ylabel('Channel')
    plt.savefig('SimWakeInDas.png',bbox_inches='tight')
