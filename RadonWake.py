
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
from math import sin,cos,tan,radians
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



def PlotRadon(ShowDataSlice,DENOISE_RADON):
    fig, (ax1, ax2) = plt.subplots(1, 2)
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
    
    if DENOISE_RADON==1:
        sinogram=(sinogram-np.mean(sinogram))/np.std(sinogram,ddof=1)
        STD=3.5*np.std(sinogram,ddof=1)
        sinogram=((sinogram>STD)|(sinogram<-STD))*sinogram

    dx, dy = 0.5 * 180.0 / max(ShowDataSlice.shape), 0.5 / sinogram.shape[0]
    #dx=0.01
    ax2.set_title("Radon transform")
    ax2.set_xlabel("Projection angle (deg)")
    ax2.set_ylabel("Projection position")
    SinSlice=slice(int(0.45*sinogram.shape[1]),int(0.55*sinogram.shape[1]),1)
    #ax2.imshow(sinogram[:,SinSlice], cmap="bwr",aspect='auto',origin='lower')
    ax2.imshow(sinogram, cmap="bwr",extent=(0.0-dx, 180.0 + dx, -dy, sinogram.shape[0] + dy),aspect='auto',origin='lower')
    plt.savefig('Radon_transform.png')
    return sinogram

def SpeedOnRadon(sinogram,resolution,channel_spacing,fs,scaling,WAVEDIRECT):
#求解最大的斜率
    
    #求解最大的斜率，平均每行的最大值索引
    #sinogram=abs(sinogram)
    deg_y_list=[]
    sino=[]
    for i in range(0,sinogram.shape[0]):
        sino.append(np.max(sinogram[i,:]))
        if WAVEDIRECT==0:  #考虑45-90
            deg_y=np.argmax(sinogram[i,0:int(sinogram.shape[1]/2)])%sinogram.shape[1]/resolution*180
            deg_y=(int(sinogram.shape[1]/4)+np.argmax(sinogram[i,int(sinogram.shape[1]/4):int(sinogram.shape[1]/2)]))%sinogram.shape[1]/resolution*180
            #print(deg_y)

        elif WAVEDIRECT==1: #考虑135-180
            deg_y=(int(sinogram.shape[1]/2)+np.argmax(sinogram[i,int(sinogram.shape[1]/2):-1]))%sinogram.shape[1]/resolution*180
        else:
            deg_y=np.argmax(sinogram[i,:])%sinogram.shape[1]/resolution*180
        deg_y_list.append(deg_y)
    df=pd.DataFrame(np.transpose(np.array([sino,deg_y_list])),columns=['sino','deg'])
    df['sino']=abs(df['sino'])
    df.sort_values(by='sino',ascending=False,inplace=True)
    deg_y_list=df['deg'][0:9]
    #print(np.mean(deg_y_list),np.median(deg_y_list))
    if tan(radians(np.mean(deg_y_list)))==0:
        s1=NaN
    else:
        s1=1/tan(radians(np.mean(deg_y_list)))*channel_spacing*fs/(scaling)
    if tan(radians(np.median(deg_y_list)))==0:
        s2=NaN
    else:
        s2=1/tan(radians(np.median(deg_y_list)))*channel_spacing*fs/(scaling)
    #print('选取Radon域上前10的最大值所对应的deg\n')
    #print('均值',s1,'中位数',s2)
    return np.mean(deg_y_list),np.median(deg_y_list),s1,s2

def SNROnRadon(sinogram):
    #[1]M. T. Rey, J. K. E. Tunaley, J. T. Folinsbee, P. A. Jahans, J. A. Dixon, and M. R. Vant, "Application Of Radon Transform Techniques To Wake Detection In Seasat-A SAR Images," IEEE Transactions on Geoscience and Remote Sensing, vol. 28, pp. 553-560, 1990.
    SNR=(np.max(sinogram)-np.mean(sinogram))/np.std(sinogram,ddof=1)
    print('SNR on Radon domain: ',SNR)
    return SNR


def ValidationSpeedOnRadon(speed,FILTER_Data,DownSampleRate,channel_spacing,MINTIME,MAXTIME,MINCHANNEL,MAXCHANNEL,Tstart,Tend,Cstart,Cend,threshold,WAVEDIRECT):
    DataCoordX, DataCoordy = FILTER_Data.shape

    TDownSample = slice(0, DataCoordX, int(1000 / DownSampleRate))
    DataDownSample = FILTER_Data[TDownSample, :]

    #画图展示的数据

    TimeWin = slice(int(MINTIME*60*DownSampleRate), max(-1,int(MAXTIME*60*DownSampleRate)), 1)
    ShowData=DataDownSample[TimeWin,slice(int(MINCHANNEL*1000/channel_spacing),int(MAXCHANNEL*1000/channel_spacing),1)]

    #对空间维度进行插值，使其接近时间的维度

    #Z-score and threshold filtering
    #ShowData=(ShowData-np.mean(ShowData))/np.std(ShowData,ddof=1)
    ShowData = stats.zscore(ShowData, axis=None)
    STD=threshold*np.std(ShowData,ddof=1)
    ShowData=((ShowData>STD)|(ShowData<-STD))*ShowData
    #Prepare data for radon transformation

    
    TimeWin = slice(int(Tstart*60*DownSampleRate), max(-1,int(Tend*60*DownSampleRate)), 1)
    Cwin = slice(int(Cstart*1000/channel_spacing),max(-1,int(Cend*1000/channel_spacing)),1)
    ShowDataSlice=(ShowData[TimeWin,Cwin])

    ShowDataSlice=np.transpose((ShowDataSlice))

    plt.figure(dpi=300)
    plt.imshow(ShowDataSlice, aspect='auto',cmap="bwr",origin='lower',vmin=-3,vmax=3)
    x=np.arange(1,ShowDataSlice.shape[1])
    k=((speed)/(DownSampleRate*channel_spacing))
    b=0
    if WAVEDIRECT==1:
        b=ShowDataSlice.shape[0]//2
    y=k*x+b
    #plt.plot(x,y,lw=3,color='g')
    plt.plot(x,k*x+1.5*b,lw=4,color='g')
    plt.text(0,0,"U =%.2f m/s"%speed,size=15)
    plt.ylim((0,ShowDataSlice.shape[0]))
    print('Validation!')
    plt.savefig('validation.png')
    return ShowDataSlice


def EnhanceResolution(ShowDataSlice,DownSampleRate,SPACE_EN,TIME_EN):
    ReDownSampleRate=DownSampleRate
    if ShowDataSlice.shape[0]>1000:
        ReDownSample=ShowDataSlice.shape[0]//1000
        TimeWin=slice(0,-1,ReDownSample)
        ShowDataSlice=ShowDataSlice[TimeWin,:]
        ReDownSampleRate=DownSampleRate/ReDownSample
    if SPACE_EN==1:
        ShowDataSlice=Dup_Spat_Dim(ShowDataSlice) #增强空间维度的细节
    if TIME_EN==1:
        ShowDataSlice=Dup_Time_Dim(ShowDataSlice) #增强时间维度的细节
    scaling=round(ShowDataSlice.shape[0]/ShowDataSlice.shape[1],4)
    channel_spacing_scaling=round((ShowDataSlice.shape[0]-1)/(ShowDataSlice.shape[1]-1),4)
    if scaling>1:
        ShowDataSlice=DASInterpolate(ShowDataSlice)  #将空间维度的尺寸补充到和时间维度一致
    else:
        scaling=1
        channel_spacing_scaling=1
    Result=ShowDataSlice
    return Result,ReDownSampleRate,channel_spacing_scaling

def CalculateKEnv(ShowDataSlice,channel_spacing,fs):
    '''
    Calculate the Sloop K of the divergent wave envelope curve
    '''
    #找出每个时刻出现振动的最大channel
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

    reg = LinearRegression()
    reg.fit(Time,Channel)
    print("The linear model is: Y = {:.5} + {:.5}X".format(reg.intercept_[0], reg.coef_[0][0]))
    K_env=reg.coef_[0][0]
    bias=reg.intercept_[0]
    plt.figure(dpi=300)
    plt.imshow(ShowDataSlice, aspect='auto',cmap="bwr",origin='lower',vmin=-3,vmax=3)
    plt.colorbar()
    plt.plot(Time,Channel,lw=2,linestyle='--',label='Envelope Curve')
    plt.plot(Time,list(bias+K_env*np.array(Time)),lw=2,label='Linear regression of Envelope Curve')
    plt.legend()
    print('Envelope!')
    plt.savefig('Envelope.png')

    return K_env,bias
