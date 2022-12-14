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
from RadonWake import PlotRadon,SpeedOnRadon,SNROnRadon,ValidationSpeedOnRadon
from scipy import stats
#%%Params ---------------------------------------
# filter or not
FILTER=1
DownSampleRate = 100    #输入的采样数据为1秒1000个点，这里设置每秒采样的点数
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

#时间差分

FILTER_Data=DataDiff(FILTER_Data,1)


#%%
#降采样
DataCoordX, DataCoordy = FILTER_Data.shape
print('DataCoordX',DataCoordX,'DataCoordy',DataCoordy)
DSR=np.arange(1,100,1)
#DSR=np.arange(250,1000,20)
#DSR=[100]
res=[]
for DownSampleRate in DSR:
# 对some_data 进行采样，因为原始数据每秒采样1000，可以降为200个点以方便人为的判断
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
    #print(ShowData)
    #Prepare data for radon transformation

    RegionSliceX=[Tstart*60*DownSampleRate,max(-1,Tend*60*DownSampleRate),max(-1,Tend*60*DownSampleRate),Tstart*60*DownSampleRate,Tstart*60*DownSampleRate]
    RegionSliceY=[int(Cstart*1000/channel_spacing),int(Cstart*1000/channel_spacing),max(-1,int(Cend*1000/channel_spacing)),max(-1,int(Cend*1000/channel_spacing)),int(Cstart*1000/channel_spacing)]
    TimeWin = slice(int(Tstart*60*DownSampleRate), max(-1,int(Tend*60*DownSampleRate)), 1)
    Cwin = slice(int(Cstart*1000/channel_spacing),max(-1,int(Cend*1000/channel_spacing)),1)
    ShowDataSlice=(ShowData[TimeWin,Cwin])

    ShowDataSlice=Dup_Spat_Dim(ShowDataSlice) #增强空间维度的细节
    ShowDataSlice=Dup_Time_Dim(ShowDataSlice) #增强时间维度的细节

    scaling=round(ShowDataSlice.shape[0]/ShowDataSlice.shape[1],4)
    channel_spacing_scaling=round((ShowDataSlice.shape[0]-1)/(ShowDataSlice.shape[1]-1),4)
    print('scaling:',scaling,ShowDataSlice.shape)
    if scaling>1:
        ShowDataSlice=DASInterpolate(ShowDataSlice)
    else:
        scaling=1
        channel_spacing_scaling=1
    print('scaling:',scaling,ShowDataSlice.shape)

    ShowDataSlice=np.transpose((ShowDataSlice))
    #%%
    PlotDAS(ShowData,ST,ET,FiberBoatMessage,MINCHANNEL,MAXCHANNEL,RegionSliceX,RegionSliceY,channel_spacing,n_channels)
    DASNR(ShowData)

    # Radon transformation and analysis
    print(ShowDataSlice.shape)
    sinogram=PlotRadon(ShowDataSlice)
    print("Radon done!")
    snr=SNROnRadon(sinogram)
    deg_mean,deg_median,s_mean,s_median=SpeedOnRadon(sinogram,max(ShowDataSlice.shape),channel_spacing,DownSampleRate,channel_spacing_scaling,WAVEDIRECT)
    res.append((DownSampleRate,snr,deg_mean,deg_median,s_mean,s_median))
    
    print(DownSampleRate)
    res.append((DownSampleRate,snr,deg_mean,deg_median,s_mean,s_median))
    


RES=pd.DataFrame(res,columns=['DownSampleRate',"snr","deg_mean","deg_median","s_mean",'s_median'])
RES.to_excel('res.xlsx')
SampleRate=list(RES['DownSampleRate'])
snr=list(RES['snr'])
Speed=list(RES['s_mean'])

fig = plt.figure(dpi=200)
ax1 = fig.add_subplot(111)

ax1.plot(SampleRate, snr)
ax1.set_ylabel('%s'%('SNR'),size=20)
plt.xlabel('Sample Rate')
ax2 = ax1.twinx()  # 设置双y轴
ax2.plot(SampleRate, Speed, 'r',)
ax2.set_ylabel('%s'%("Gradient K"),size=20,color='r')
plt.savefig('SNR.png')
#%%

#波线斜率
RES.dropna(inplace=True)
RES.sort_values(by='snr',ascending=False,inplace=True)
#
SNR_threshold=np.mean(RES['snr'])+0.75*(np.max(RES['snr'])-np.mean(RES['snr']))

RES1=RES[RES['snr']>SNR_threshold]
s_mean = RES1['s_mean']

print('波线斜率',np.mean(s_mean))

#在ShowDataslice上补充波线斜率直线
speed=30
DownSampleRate = 1000
ValidationSpeedOnRadon(speed,FILTER_Data,DownSampleRate,channel_spacing,MINTIME,MAXTIME,MINCHANNEL,MAXCHANNEL,Tstart,Tend,Cstart,Cend,threshold)