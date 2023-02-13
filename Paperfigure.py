from ast import And
from cProfile import label
from cmath import exp
from email import header
from operator import and_
from pyexpat import features
from tracemalloc import start
from turtle import color
import preprocessing
import tkinter
import pandas as pd
import numpy as np
import os
import time
import matplotlib as mpl
from scipy.misc import central_diff_weights
import statsmodels.api as sm
import statsmodels.formula.api as smf
import pmdarima as pm
from numpy import array, sign, zeros
from scipy.interpolate import interp1d
import scipy.signal
import math
from math import sin,cos,tan,radians,sqrt
from sklearn import preprocessing
import matplotlib.pyplot as plt
import scipy.signal as signal
from scipy.fftpack import fft,ifft
from sklearn import preprocessing
from sklearn.decomposition import FastICA
import datetime
from datetime import timedelta
from DASFileRead import DasFileRead
from tdms_reader_1 import *
from AISData import AISData,AnchorShip
from DASFilter import bandpass_f,WeightMatrix
from numpy import NaN
import skimage 
from skimage.transform import radon
from sklearn.linear_model import LinearRegression
from scipy import stats
from DASFilter import bandpass_f,DataDiff,WeightMatrix,PlotDAS,DASNR,DASInterpolate,Dup_Spat_Dim,Dup_Time_Dim
from WaveVAngel import PlotSimuWaveInDas,FroudeNum,WavePattern,ROTATE,move,cross,WavePattern1



def PlotRadonInPaper(ShowDataSlice):
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    plt.subplots_adjust(wspace=0.5)
    ax1.set_title("Original")
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Channel')
    print('ShowDataSlice size',ShowDataSlice.shape)
    ax1.imshow(ShowDataSlice, aspect='auto',cmap="bwr",origin='lower',vmin=-3,vmax=3)
    theta = np.linspace(0., 180., max(ShowDataSlice.shape), endpoint=False)
    #theta=np.linspace(0.0,180.0,12,endpoint=False)
    sinogram = radon(ShowDataSlice, theta=theta)
    print('Radon size',sinogram.shape)
    
    dx, dy = 0.5 * 180.0 / max(ShowDataSlice.shape), 0.5 / sinogram.shape[0]
    #dx=0.01
    ax2.set_title("Radon transform")
    ax2.set_xlabel("Projection angle (deg)")
    ax2.set_ylabel("Projection position")
    SinSlice=slice(int(0.45*sinogram.shape[1]),int(0.55*sinogram.shape[1]),1)
    #ax2.imshow(sinogram[:,SinSlice], cmap="bwr",aspect='auto',origin='lower')
    ax2.imshow(sinogram, cmap="bwr",extent=(0.0-dx, 180.0 + dx, -dy, sinogram.shape[0] + dy),aspect='auto',origin='lower')


    sinogram=(sinogram-np.mean(sinogram))/np.std(sinogram,ddof=1)
    STD=3.5*np.std(sinogram,ddof=1)
    sinogram=((sinogram>STD)|(sinogram<-STD))*sinogram

    dx, dy = 0.5 * 180.0 / max(ShowDataSlice.shape), 0.5 / sinogram.shape[0]
    #dx=0.01
    ax3.set_title("Radon transform\n(denoised)")
    ax3.set_xlabel("Projection angle (deg)")
    ax3.set_ylabel("Projection position")
    SinSlice=slice(int(0.45*sinogram.shape[1]),int(0.55*sinogram.shape[1]),1)
    #ax2.imshow(sinogram[:,SinSlice], cmap="bwr",aspect='auto',origin='lower')
    ax3.imshow(sinogram, cmap="bwr",extent=(0.0-dx, 180.0 + dx, -dy, sinogram.shape[0] + dy),aspect='auto',origin='lower')
    plt.savefig('Radon_transform.png' ,bbox_inches='tight')
    return sinogram


def PlotK_KenvLine(ShowDataSlice,speed,DownSampleRate,channel_spacing,WAVEDIRECT,K_env,bias):
    '''
    depict the divergent wave and the envolope curve based on the estimated speed in a data slice
    '''
    plt.figure(dpi=300)
    plt.imshow(ShowDataSlice, aspect='auto',cmap="bwr",origin='lower',vmin=-3,vmax=3)
    x=np.arange(1,ShowDataSlice.shape[1])
    k=((speed)/(DownSampleRate*channel_spacing))
    b=0
    if WAVEDIRECT==1:
        b=ShowDataSlice.shape[0]//2
    y=k*x+b
    #plt.plot(x,y,lw=3,color='g')
 
    Time=[]
    Channel=[]
    for i in range(0,ShowDataSlice.shape[1]):
        Time.append(i)
        for j in range(ShowDataSlice.shape[0]-1,0,-1):
            if ShowDataSlice[j,i]!=0:
                Channel.append(j)
                break
    Time=np.array(Time).reshape(-1,1)
    Channel=np.array(Channel).reshape(-1,1)

    plt.figure(dpi=300)
    plt.plot(x,k*x+1.5*b,lw=2,color='lime',label='Calculated divergent wave crest')
    #plt.text(0,0,"U =%.2f m/s"%speed,size=15)
    plt.ylim((0,ShowDataSlice.shape[0]))

    plt.imshow(ShowDataSlice, aspect='auto',cmap="bwr",origin='lower',vmin=-3,vmax=3)
    plt.plot(Time,Channel,lw=2,linestyle='--',color='b',label='Envelope Curve')
    plt.plot(Time,list(bias+K_env*np.array(Time)),lw=2,color='gold',label='Linear regression of Envelope Curve')
    plt.legend()
    print('Envelope and wave crest!')
    plt.savefig('Envelope_Crest.png' ,bbox_inches='tight')


def PlotSimulInDAS(DownSampleRate,v,h,angle,A,ShowData,ST,ET,MINCHANNEL,channel_spacing,WLen_Scale,Wbias,Tbias,UpperBound,LowerBound):
    '''
    Depict the simulated ship wake in the DAS figure to analyze the difference between the simulated  and measured divergent wave
    '''
    T=30
    N=FroudeNum(v,h)
    print('Froude Number is ',N)



    delta_T=0.1
    plt.figure(dpi=300,figsize=(10,10))
    for a in range(1,A+1):
        crossp1=[]
        crossp2=[]
        t_start=0
        t_start1=0
        for t in np.arange(0,T,delta_T):
            dist=v*t #移动的距离 dist=v*t,除以矫正系数
            X,Y,Y1,ALPHA=WavePattern1(a,N,0.2,1.5,0.01)
            alpha=max(ALPHA)
            Attenuation0 = sin(radians(angle - (alpha)))**2   #乘以一定系数
            Attenuation1 = sin(radians(180-angle - (alpha)))**2
            RX,RY=ROTATE(angle,X,Y)
            RX1,RY1=ROTATE(angle,X,Y1)
            RMX,RMY=move(RX,RY,angle,dist)
            RMX1,RMY1=move(RX1,RY1,angle,dist)
            cp=cross(RMX,RMY)
            #Apply the upper and the lower bound of the wave crest to fit the DAS data
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
        x1=x1*delta_T+t_start+Tbias
        x2=np.arange(0,len(crossp2),1)
        x2=x2*delta_T+t_start1+Tbias
        
        #Scale the wave length because the WavePattern function returns the dimensionless wake pattern
        x1=WLen_Scale*x1
        x2=WLen_Scale*x2
        crossp1=WLen_Scale*(np.array(crossp1)+Wbias)
        crossp2=WLen_Scale*(np.array(crossp2)+Wbias)

        #Channel sampling and sampling of the simulated wave crest
        crossp1=crossp1/channel_spacing
        crossp2=crossp2/channel_spacing
        x1=x1*DownSampleRate
        x2=x2*DownSampleRate
        plt.plot(x1,crossp1,lw=2.5,linestyle='--',alpha=Attenuation0,color='lime')
        plt.plot(x2,crossp2,lw=2.5,linestyle='--',alpha=Attenuation1,color='lime')

    print(np.transpose(ShowData).shape)
    plt.imshow(np.transpose(ShowData), cmap="bwr", aspect='auto',origin='lower',vmin=-3,vmax=3) # cmap=''bwr,, 
    #计算坐标轴的刻度大小,合计10个时间戳(计算有问题，需要考虑数据的实际距离以及截断)
    #根据MINCHANNEL和MAXCHANNEL重置Y轴
    plt.colorbar()
    TimeTicks=10
    TI=(ET-ST)/TimeTicks
    xlabel=np.linspace(0,ShowData.shape[0],TimeTicks+1)
    plt.xticks(xlabel,pd.date_range(ST.strftime("%Y%m%d %H%M"),ET.strftime("%Y%m%d %H%M"),freq=TI).strftime('%H:%M:%S'),rotation = 60,size=15)
    ylabel=np.linspace(0,ShowData.shape[1],10)
    plt.xlabel("Time",fontsize=15)
    plt.yticks(ylabel,np.round((ylabel)*channel_spacing/1000+MINCHANNEL,2),size=15)

    plt.ylabel("Distance(Km)",fontsize=15)
    plt.savefig('PlotSimulInDAS.png',bbox_inches='tight')