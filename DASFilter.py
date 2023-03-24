
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
from math import sin,cos
#mpl.use('tKAgg')
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
from scipy import interpolate

def bandpass_f(some_data,sample_rate,cutoff_f_low,cutoff_f_high,cutoff_order):
    #带通滤波
    data=[]
    cutoff_f0=cutoff_f_low/sample_rate*2 #截止频率10hz
    cutoff_f1=cutoff_f_high/sample_rate*2 
    #cutoff_order=cutoff_order   #配置滤波器的阶数
    if cutoff_f0==0 and cutoff_f1!=0:
        b,a = signal.butter(cutoff_order, cutoff_f1, 'low')
    elif cutoff_f0!=0 and cutoff_f1==0:
        b,a = signal.butter(cutoff_order, cutoff_f0, 'high')
    else:
        b,a = signal.butter(cutoff_order, [cutoff_f0,cutoff_f1], 'bandpass') 
    #分别对所有的通道滤波，滤波后的数据保存在some_data里
    for channel in range(0,some_data.shape[1]):
        data.append(signal.filtfilt(b,a, some_data[:,channel]))
        if (some_data.shape[1]-channel)%100==0:
            print('The rest channels for filtering:{}'.format(some_data.shape[1]-channel))
    return np.array(data).transpose()

def DataDiff(Data,order):
    for i in range(1,order):
        Data=np.diff(Data,axis=0)
    return Data


def WeightMatrix(K,size):
    Weight_matrix=np.zeros((size,size))
    X=np.arange(0,size,1)
    b=size/2-K*size/2
    Y=K*X+b
    Y=np.round(Y,0)
    Cdim=[]  #Channel dimension
    Tdim=[]  #Time dimension
    for x in X:
        if int(Y[x]) <size and int(Y[x])>=0:
            Weight_matrix[size-1-int(Y[x])][int(x)]=1
            Cdim.append(size-1-int(Y[x]))
            Tdim.append(int(x))

    return Weight_matrix,Cdim,Tdim


def PlotDAS(ShowData,ST,ET,FiberBoatMessage,MINCHANNEL,MAXCHANNEL,RegionSliceX,RegionSliceY,channel_spacing,n_channels,PLOTANCHOR,PLOTREGION):
#plot figure for ShowData, and mark the passing ship according to the message from FiberBoatMessage. Besides, the data slice used for the next radon transfromation is also mark in this figure
    fig1=plt.figure(dpi=400,figsize=(10,6))    
    ax1 = fig1.add_subplot(1,1,1) 
    plt.imshow(np.transpose(ShowData), cmap="bwr", aspect='auto',origin='lower',vmin=-3,vmax=3) # cmap=''bwr,, 
    #计算坐标轴的刻度大小,合计10个时间戳(计算有问题，需要考虑数据的实际距离以及截断)
    #根据MINCHANNEL和MAXCHANNEL重置Y轴
    plt.colorbar()
    TimeTicks=3
    TI=(ET-ST)/TimeTicks
    xlabel=np.linspace(0,ShowData.shape[0]-1,TimeTicks+1)
    plt.xticks(xlabel,pd.date_range(ST.strftime("%Y%m%d %H%M%S"),ET.strftime("%Y%m%d %H%M%S"),periods=TimeTicks+1).strftime('%H:%M:%S'),rotation = 0,size=15)
    ylabel=np.linspace(0,ShowData.shape[1],5)
    plt.xlabel("Time",fontsize=15)
    plt.yticks(ylabel,np.round((ylabel)*channel_spacing/1000+MINCHANNEL,2),size=15)

    plt.ylabel("Distance(Km)",fontsize=15)
    #船只过光纤轨迹标定
    
    if PLOTANCHOR==1:
        deltaCTUpper,deltaCTdown,RegionDistUpper,RegionDistdown=AnchorShip(FiberBoatMessage,MINCHANNEL,MAXCHANNEL,n_channels,channel_spacing,ST,ET,ShowData.shape[0])
        print(RegionDistUpper,RegionDistdown)
        for i in range(0,len(deltaCTdown)):
            plt.plot([deltaCTdown[i],deltaCTUpper[i],deltaCTUpper[i],deltaCTdown[i],deltaCTdown[i]],[RegionDistdown[i],RegionDistdown[i],   RegionDistUpper[i],RegionDistUpper[i],RegionDistdown[i]],linewidth=1)
    if PLOTREGION==1:
        plt.plot(RegionSliceX,RegionSliceY,linewidth=2.5)
    plt.savefig("Space-time diagram.png",bbox_inches = 'tight')
    plt.savefig("/home/huangwj/DAS/BoatTrajectory/Paperfig/Space-time diagram2.pdf",bbox_inches = 'tight')
 
    print("Space-time diagram.png")

#def DASNR():
def DASNR(ShowData):
    #[1]M. T. Rey, J. K. E. Tunaley, J. T. Folinsbee, P. A. Jahans, J. A. Dixon, and M. R. Vant, "Application Of Radon Transform Techniques To Wake Detection In Seasat-A SAR Images," IEEE Transactions on Geoscience and Remote Sensing, vol. 28, pp. 553-560, 1990.
    SNR=(np.max(ShowData)-np.mean(ShowData))/np.std(ShowData,ddof=1)
    print('SNR on image domain: ',SNR)
    return SNR


def DASInterpolate(ShowData):
    '''
    空间维度插值，使其等效于时间维度
    '''
    SD=ShowData.shape[1]
    TD=ShowData.shape[0]
    Y=[]
    for t in range(0,TD):
        f = interpolate.interp1d(np.arange(0,SD,1), ShowData[t,:], kind='linear')
        x=list(np.linspace(0,SD-1,TD,endpoint=True))
        y=list(np.round(f(x),4))
        Y.append(y)
    return np.array(Y)

def Dup_Spat_Dim(ShowData):
    if ShowData.shape[1]<1000:
        loop=1000//ShowData.shape[1]
        temp=ShowData    
        for i in range(0,loop-1):
            temp=np.hstack((temp,ShowData))
        ShowData=temp
    return ShowData


def Dup_Time_Dim(ShowData):
    if ShowData.shape[0]<1000:
        loop=1000//ShowData.shape[0]
        temp=ShowData    
        for i in range(0,loop-1):
            temp=np.vstack((temp,ShowData))
        ShowData=temp
    return ShowData