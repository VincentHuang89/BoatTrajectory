#读取DAS数据保存在/home/huangwj/DAS/BoatTrajectory/Data
#%%
from ast import And
from cProfile import label
from cmath import exp
from email import header
from operator import and_
from pyexpat import features
from tracemalloc import start
from turtle import color

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
from sklearn.decomposition import FastICA
import datetime
from datetime import timedelta
from DASFileRead import DasFileRead
from tdms_reader_1 import *
from AISData import AISData,AnchorShip
from DASFilter import bandpass_f,WeightMatrix,PlotDAS
import skimage 
from skimage.transform import radon
from RadonWake import PlotRadon,SpeedOnRadon,SNROnRadon

#%%Params ---------------------------------------
# filter or not
FILTER=1
DownSampleRate = 100    #输入的采样数据为1秒1000个点，这里设置每秒采样的点数
# Showdata params
MINTIME=0
MAXTIME=-1
MINCHANNEL=0
MAXCHANNEL=1500
# Z-score and threshold filter
threshold=2.5
#radon transfromation params
Tstart =2 #2.5  #minute
Tend =3  #5.5   #minute
Cstart = 2.5 #1800  #Km
Cend= 3.6 #2500    #Km

#---------------------------------------------------

start = time.time()
mpl.rcParams['font.sans-serif'] = ['SimHei']  # 显示中文
mpl.rcParams['axes.unicode_minus'] = False  # 显示负号
# from func import *
# 设定工作路径
path = '/home/huangwj/DAS/BoatTrajectory/'
os.chdir(path)
print('当前工作路径是：' + os.getcwd())  # 显示当前路径


#选择数据路径
'''
对应的船只MMSI和轨迹时间段（UTC+8）分别为
413702130 |  2022-08-08: 09:00-09:15  需要扩充到9点30分
352929000 |  2022-08-08: 12:37-12:47
636017493 |  2022-08-08: 18:05-18:20
413478650 |  2022-08-11: 14:22-15:00 
413480930 |  2022-08-11: 23:30- 2022-08-12: 00:05
tdms文件的时间是以UTC+0来命名，两者存在区别
'''
DataPath='/home/huangwj/DAS/BoatTrajectory/DataforAIS'
ST_UTC8=datetime.datetime.strptime("27/07/22 8:00", "%d/%m/%y %H:%M")
ET_UTC8=datetime.datetime.strptime("27/07/22 8:59", "%d/%m/%y %H:%M")
ST_UTC0=ST_UTC8-timedelta(hours=8)
ET_UTC0=ET_UTC8-timedelta(hours=8)
FileSet,times=DasFileRead(ST_UTC0,ET_UTC0,DataPath)
ST=times[0]
ET=times[-1]+timedelta(minutes=1)
ST=ST+timedelta(hours=8)
ET=ET+timedelta(hours=8)
print('Minutes of the DAS data is ',len(FileSet))

#%%读取文件
for i in range(0,len(FileSet)):
    FilePath= os.path.join(DataPath,FileSet[i])

    if i==0:
        tdms = TdmsReader(FilePath)
        props = tdms.get_properties()
        zero_offset = props.get('Zero Offset (m)')
        channel_spacing = props.get('SpatialResolution[m]') * props.get('Fibre Length Multiplier')
        n_channels = tdms.fileinfo['n_channels']
        depth = zero_offset + np.arange(n_channels) * channel_spacing
        fs = props.get('SamplingFrequency[Hz]')
        print('Channel_spacing is:{:.4f}m'.format(channel_spacing))
        print('Number of channels in file: {0}'.format(n_channels))
        print('Time samples in file: {0}'.format(tdms.channel_length))
        print('Sampling frequency (Hz): {0}'.format(fs))
        first_channel = 0  
        last_channel = n_channels - 1  
        first_time_sample = 00 
        last_time_sample = tdms.channel_length - 1 
        some_data = tdms.get_data(first_channel, last_channel, first_time_sample, last_time_sample)   
        Data = some_data
    else:
        tdms = TdmsReader(FilePath)
        some_data1 = tdms.get_data(first_channel, last_channel, first_time_sample, last_time_sample)
        Data = np.concatenate((Data, some_data1), axis=0)
    print(FilePath)

#%%
#读取时间区间的船只数据
PosFile='pos_zf_gsd_1657814400_1660579199_742.csv'
StaticFile='static_zf_gsd_1657814400_1660579199_742.csv'
FiberBoatMessage=AISData(PosFile,StaticFile,ST,ET)

#%%滤波

if FILTER==1:
    FILTER_Data = bandpass_f(Data, fs, 0,0.5,4) 
    print('Filter!')
else:
    FILTER_Data=Data

#%%
#降采样
DataCoordX, DataCoordy = FILTER_Data.shape
print('DataCoordX',DataCoordX,'DataCoordy',DataCoordy)
# 对some_data 进行采样，因为原始数据每秒采样1000，可以降为200个点以方便人为的判断
TDownSample = slice(0, DataCoordX, int(1000 / DownSampleRate))
DataDownSample = FILTER_Data[TDownSample, :]


#画图展示的数据

TimeWin = slice(int(MINTIME*60*DownSampleRate), max(-1,int(MAXTIME*60*DownSampleRate)), 1)
ShowData=DataDownSample[TimeWin,slice(MINCHANNEL,MAXCHANNEL,1)]

#Z-score and threshold filtering
ShowData=(ShowData-np.mean(ShowData))/np.std(ShowData,ddof=1)
STD=threshold*np.std(ShowData,ddof=1)
ShowData=((ShowData>STD)|(ShowData<-STD))*ShowData

#Prepare data for radon transformation

RegionSliceX=[Tstart*60*DownSampleRate,Tend*60*DownSampleRate,Tend*60*DownSampleRate,Tstart*60*DownSampleRate,Tstart*60*DownSampleRate]
RegionSliceY=[int(Cstart*1000/channel_spacing),int(Cstart*1000/channel_spacing),int(Cend*1000/channel_spacing),int(Cend*1000/channel_spacing),int(Cstart*1000/channel_spacing)]
TimeWin = slice(int(Tstart*60*DownSampleRate), int(Tend*60*DownSampleRate), 1)
Cwin = slice(int(Cstart*1000/channel_spacing),int(Cend*1000/channel_spacing),1)
ShowDataSlice=(ShowData[TimeWin,Cwin])
ShowDataSlice=np.transpose((ShowDataSlice))
#%%
PlotDAS(ShowData,ST,ET,FiberBoatMessage,MINCHANNEL,RegionSliceX,RegionSliceY,channel_spacing,n_channels)
# Radon transformation and analysis
sinogram=PlotRadon(ShowDataSlice)
SNROnRadon(sinogram)
SpeedOnRadon(sinogram,max(ShowDataSlice.shape),channel_spacing,DownSampleRate)
#%%
