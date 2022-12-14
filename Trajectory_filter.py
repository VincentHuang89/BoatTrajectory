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
from DASFilter import bandpass_f,DataDiff,WeightMatrix,PlotDAS,DASNR,DASInterpolate,Dup_Spat_Dim,Dup_Time_Dim
import skimage 
from skimage.transform import radon
from RadonWake import PlotRadon,SpeedOnRadon,SNROnRadon,ValidationSpeedOnRadon,EnhanceResolution
from scipy import stats
#%%Params ---------------------------------------
# filter or not
FILTER=1
DownSampleRate = 50   #输入的采样数据为1秒1000个点，这里设置每秒采样的点数
# Showdata params
MINTIME=0
MAXTIME=-1
MINCHANNEL=0
MAXCHANNEL=14   #Km
#波线方向（左下到右上：0 （图像域的上半部分，deg：0-90），左上到右下：1，deg：90-180，展示所有：2）
WAVEDIRECT=0

# Z-score and threshold filter
threshold=2
#radon transfromation params
Tstart =1.8 #2.8  #1.8   #minute
Tend =2# 3 #2     #minute
Cstart =11.15#10 #11.15 #1800  #Km
Cend= 11.3#11 #11.3 #2500    #Km

Cstart=max(Cstart-MINCHANNEL,0)
Cend=min(Cend-MINCHANNEL,MAXCHANNEL)
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
ST_UTC8=datetime.datetime.strptime("24/07/22 21:08", "%d/%m/%y %H:%M")
ET_UTC8=datetime.datetime.strptime("24/07/22 21:11", "%d/%m/%y %H:%M")
ST_UTC0=ST_UTC8-timedelta(hours=8)
ET_UTC0=ET_UTC8-timedelta(hours=8)
FileSet,times=DasFileRead(ST_UTC0,ET_UTC0,DataPath)
ST=times[0]
ET=times[-1]+timedelta(minutes=1)
ST=ST+timedelta(hours=8)
ET=ET+timedelta(hours=8)
print('Minutes of the DAS data is ',len(FileSet))

#%%Read DAS data files
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
#Read AIS data files
PosFile='pos_zf_gsd_1657814400_1660579199_742.csv'
StaticFile='static_zf_gsd_1657814400_1660579199_742.csv'
FiberBoatMessage=AISData(PosFile,StaticFile,ST,ET)

#%%Filter the raw data

if FILTER==1:
    FILTER_Data = bandpass_f(Data, fs, 0,0.5,4) 
    print('Filter!')
else:
    FILTER_Data=Data

#Get the differential data in the time dimension
FILTER_Data=DataDiff(FILTER_Data,1)


#%%
#Downsample FILTER_Data
DataCoordX, DataCoordy = FILTER_Data.shape
print('DataCoordX',DataCoordX,'DataCoordy',DataCoordy)
TDownSample = slice(0, DataCoordX, int(1000 / DownSampleRate))
DataDownSample = FILTER_Data[TDownSample, :]

#Obtain the data according to the Time and Channel Parameters from the downsampled data, to show the DAS image
TimeWin = slice(int(MINTIME*60*DownSampleRate), max(-1,int(MAXTIME*60*DownSampleRate)), 1)
ShowData=DataDownSample[TimeWin,slice(int(MINCHANNEL*1000/channel_spacing),int(MAXCHANNEL*1000/channel_spacing),1)]

#Z-score and threshold filtering
ShowData = stats.zscore(ShowData, axis=None)
STD=threshold*np.std(ShowData,ddof=1)
ShowData=((ShowData>STD)|(ShowData<-STD))*ShowData

#plot ShowData and anchor the data sliced region
RegionSliceX=[Tstart*60*DownSampleRate,max(-1,Tend*60*DownSampleRate),max(-1,Tend*60*DownSampleRate),Tstart*60*DownSampleRate,Tstart*60*DownSampleRate]
RegionSliceY=[int(Cstart*1000/channel_spacing),int(Cstart*1000/channel_spacing),max(-1,int(Cend*1000/channel_spacing)),max(-1,int(Cend*1000/channel_spacing)),int(Cstart*1000/channel_spacing)]
PlotDAS(ShowData,ST,ET,FiberBoatMessage,MINCHANNEL,MAXCHANNEL,RegionSliceX,RegionSliceY,channel_spacing,n_channels)
DASNR(ShowData)


#For Radon transformation, Slice partial data from ShowData-----------------------------
TimeWin = slice(int(Tstart*60*DownSampleRate), max(-1,int(Tend*60*DownSampleRate)), 1)
Cwin = slice(int(Cstart*1000/channel_spacing),max(-1,int(Cend*1000/channel_spacing)),1)
ShowDataSlice=(ShowData[TimeWin,Cwin])

#readjust the time dimension of the ShowDataSlice, to avoid the time-consuming calcualtion in the radon transformation caused by the too long time dimension
ShowDataSlice,ReDownSampleRate,channel_spacing_scaling=EnhanceResolution(ShowDataSlice,DownSampleRate)
print('ReDownSampleRate',ReDownSampleRate)
ShowDataSlice=np.transpose((ShowDataSlice))

# Radon transformation and analysis
sinogram=PlotRadon(ShowDataSlice)
print("Radon done!")
snr=SNROnRadon(sinogram)
deg_mean,deg_median,s_mean,s_median=SpeedOnRadon(sinogram,max(ShowDataSlice.shape),channel_spacing,ReDownSampleRate,channel_spacing_scaling,WAVEDIRECT)
res=[]
res.append((ReDownSampleRate,snr,deg_mean,deg_median,s_mean,s_median))
RES=pd.DataFrame(res,columns=['DownSampleRate',"snr","deg_mean","deg_median","s_mean",'s_median'])
RES.to_excel('res.xlsx')
RES.dropna(inplace=True)
s_mean = RES['s_mean']
speed=np.mean(s_mean)
print('Estimated ship speed: ',speed)

#To validate the accuracy of the estimated speed, plot the line according to the estimated speed in the ShowDataSlice image.
ValidationSpeedOnRadon(speed,FILTER_Data,ReDownSampleRate,channel_spacing,MINTIME,MAXTIME,MINCHANNEL,MAXCHANNEL,Tstart,Tend,Cstart,Cend,threshold)