
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
import pywt
from math import sqrt,log
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


def DWT(Signal,wavelet,level):
    coeffs = pywt.wavedec(Signal, wavelet, level=level)
    thresh  = np.median(np.abs(Signal))/0.6745*sqrt(2*log(len(Signal)))
    coeffs[1:] = (pywt.threshold(i, value=thresh, mode='soft') for i in coeffs[1:])
    return pywt.waverec(coeffs, wavelet)

def DWT_all_channels(some_data,wavelet,level):
    data=[]
    for channel in range(0,some_data.shape[1]):
        data.append(DWT(some_data[:,channel],wavelet,level))
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


def PlotDAS(ShowData,ST,ET,FiberBoatMessage,MINCHANNEL,MAXCHANNEL,REGION,channel_spacing,n_channels,zero_offset,PLOTANCHOR,PLOTREGION,SHIP):
#plot figure for ShowData, and mark the passing ship according to the message from FiberBoatMessage. Besides, the data slice used for the next radon transfromation is also mark in this figure
    ShipWave_top_pos = REGION['Wave_peak']
    fig1=plt.figure(dpi=400,figsize=(15,6))    
    #ax1 = fig1.add_subplot(1,1,1) 
    plt.imshow(np.transpose(ShowData), cmap="bwr", aspect='auto',origin='lower',vmin=-3,vmax=3,extent=[0, ShowData.shape[0], 0, ShowData.shape[1]]) # cmap=''bwr,
    plt.colorbar()

    #plt.scatter(ShipWave_top_pos[0], ShipWave_top_pos[1], color='yellow', marker='o')
    #计算坐标轴的刻度大小,合计10个时间戳(计算有问题，需要考虑数据的实际距离以及截断)
    #根据MINCHANNEL和MAXCHANNEL重置Y轴
    TimeTicks = 7
    xlabel = np.linspace(0, ShowData.shape[0] - 1, TimeTicks)
    xtick_labels = pd.date_range(ST.strftime("%Y%m%d %H%M%S"), ET.strftime("%Y%m%d %H%M%S"), periods=TimeTicks).strftime('%H:%M:%S')
    TimeTicks=7
    xlabel=np.linspace(0,ShowData.shape[0]-1,TimeTicks)
    plt.xticks(xlabel,xtick_labels,rotation = 0,size=15)
    ylabel=np.linspace(0,ShowData.shape[1],5)
    plt.xlabel("Time",fontsize=15)
    plt.yticks(ylabel,np.round((ylabel)*channel_spacing/1000+MINCHANNEL,2),size=15)
    plt.ylabel("Distance(Km)",fontsize=15)

    #船只过光纤轨迹标定
    
    if PLOTANCHOR==1:
        deltaCTUpper,deltaCTdown,RegionDistUpper,RegionDistdown,_=AnchorShip(FiberBoatMessage,MINCHANNEL,MAXCHANNEL,n_channels,channel_spacing,zero_offset,ST,ET,ShowData.shape[0],ShowData.shape[1])
        
        for i in range(0,len(deltaCTdown)):
            plt.plot([deltaCTdown[i],deltaCTUpper[i],deltaCTUpper[i],deltaCTdown[i],deltaCTdown[i]],[RegionDistdown[i],RegionDistdown[i],   RegionDistUpper[i],RegionDistUpper[i],RegionDistdown[i]],linewidth=1)
    
    if PLOTREGION==1:
        Region_all_pos =REGION['K1']
        RegionSliceX = [Region_all_pos[0][0],Region_all_pos[1][0],Region_all_pos[2][0],Region_all_pos[3][0],Region_all_pos[0][0]]
        RegionSliceY = [Region_all_pos[0][1],Region_all_pos[1][1],Region_all_pos[2][1],Region_all_pos[3][1],Region_all_pos[0][1]]        
        plt.plot(RegionSliceX,RegionSliceY,linewidth=2.5)

        Region_all_pos =REGION['K2']
        RegionSliceX = [Region_all_pos[0][0],Region_all_pos[1][0],Region_all_pos[2][0],Region_all_pos[3][0],Region_all_pos[0][0]]
        RegionSliceY = [Region_all_pos[0][1],Region_all_pos[1][1],Region_all_pos[2][1],Region_all_pos[3][1],Region_all_pos[0][1]]        
        plt.plot(RegionSliceX,RegionSliceY,linewidth=2.5)

        Region_all_pos =REGION['Ke']
        RegionSliceX = [Region_all_pos[0][0],Region_all_pos[1][0],Region_all_pos[2][0],Region_all_pos[3][0],Region_all_pos[0][0]]
        RegionSliceY = [Region_all_pos[0][1],Region_all_pos[1][1],Region_all_pos[2][1],Region_all_pos[3][1],Region_all_pos[0][1]]        
        plt.plot(RegionSliceX,RegionSliceY,linewidth=2.5)

    plt.savefig(f"Space-time diagram_{SHIP}.svg",bbox_inches = 'tight')
    plt.savefig(f"/home/huangwj/DAS/BoatTrajectory/Paperfig/Space-time diagram_{SHIP}.pdf",bbox_inches = 'tight')
 
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

def Cal_pos_t(Ship_relocate_time,DownSampleRate):
    '''
    已知船行波尖峰时刻到图片开始时间之间的差值，计算尖峰的位置
    '''
    pos_t = int(DownSampleRate*Ship_relocate_time.total_seconds())
    return pos_t

def Cal_pos_s(POS_in_space,channel_spacing):
    '''
    已知船行波尖峰位置与图片的最低位置之间的差值，计算减分的位置
    '''
    pos_s = int((POS_in_space)*1000/channel_spacing)
    return pos_s



def Pos_AnchorRegion(Ship_relocate_time,Ship_relocate_space,DownSampleRate,channel_spacing):
    '''
    返回船行波尖峰在ShowData矩阵中的位置
    '''
    pos_t = Cal_pos_t(Ship_relocate_time,DownSampleRate)
    pos_s = Cal_pos_s(Ship_relocate_space,channel_spacing)
    ShipWave_top_pos = (pos_t,pos_s)
    return ShipWave_top_pos


def Cal_region_all_pos(ShipWave_top_pos,Bias_region_pos,Len_region_pos):
    '''
    返回切片区域在showdata矩阵中的四个顶点坐标
    '''
    left_up_point = (ShipWave_top_pos[0]+Bias_region_pos[0],ShipWave_top_pos[1]+Bias_region_pos[1])
    right_up_point = (ShipWave_top_pos[0]+Bias_region_pos[0]+Len_region_pos[0],ShipWave_top_pos[1]+Bias_region_pos[1])
    left_down_point = (ShipWave_top_pos[0]+Bias_region_pos[0],ShipWave_top_pos[1]+Bias_region_pos[1]-Len_region_pos[1])
    right_down_point = (ShipWave_top_pos[0]+Bias_region_pos[0]+Len_region_pos[0],ShipWave_top_pos[1]+Bias_region_pos[1]-Len_region_pos[1])
    return [left_up_point,right_up_point,right_down_point,left_down_point]





def Params_ShowData_slice(Ship_relocate,ST,MinChannel,Bias,Len_region,DownSampleRate,channel_spacing):
    Ship_relocate_time = datetime.datetime.strptime(Ship_relocate[0],"%d/%m/%y %H:%M:%S")-ST
    Ship_relocate_space = Ship_relocate[1]-MinChannel
    ShipWave_top_pos = Pos_AnchorRegion(Ship_relocate_time,Ship_relocate_space,DownSampleRate,channel_spacing)
    Bias_region_pos = (Cal_pos_t(Bias[0]*timedelta(minutes=1),DownSampleRate), Cal_pos_s(Bias[1],channel_spacing))
    Len_region_pos = (Cal_pos_t(Len_region[0]*timedelta(minutes=1),DownSampleRate),Cal_pos_s(Len_region[1],channel_spacing))
    Region_all_pos = Cal_region_all_pos(ShipWave_top_pos,Bias_region_pos,Len_region_pos)
    return Region_all_pos,ShipWave_top_pos

def Params_ShowData_3slice(Ship_relocate,ST,MinChannel,Bias,Len_region,DownSampleRate,channel_spacing,SHIP):
    Bias_K1 = Bias['k1'][SHIP]
    Len_region_K1 = Len_region['k1'][SHIP]
    Ship_relocate = Ship_relocate[SHIP]
    Region_all_pos_K1,ShipWave_top_pos = Params_ShowData_slice(Ship_relocate,ST,MinChannel,Bias_K1,Len_region_K1,DownSampleRate,channel_spacing)
    Bias_K2 = Bias['k2'][SHIP]
    Len_region_K2 = Len_region['k2'][SHIP]
    Region_all_pos_K2,_ = Params_ShowData_slice(Ship_relocate,ST,MinChannel,Bias_K2,Len_region_K2,DownSampleRate,channel_spacing)
    Bias_Ke = Bias['ke'][SHIP]
    Len_region_Ke = Len_region['ke'][SHIP]
    Region_all_pos_Ke,_ = Params_ShowData_slice(Ship_relocate,ST,MinChannel,Bias_Ke,Len_region_Ke,DownSampleRate,channel_spacing)
    #打包船行波切片数据的设置
    REGION = {'Wave_peak': ShipWave_top_pos,'K1':Region_all_pos_K1,'K2':Region_all_pos_K2,'Ke':Region_all_pos_Ke} 
    return REGION
