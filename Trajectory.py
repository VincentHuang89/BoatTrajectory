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
from AISData import AISData
#%%
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
ST_UTC8=datetime.datetime.strptime("24/07/22 10:05", "%d/%m/%y %H:%M")
ET_UTC8=datetime.datetime.strptime("24/07/22 10:20", "%d/%m/%y %H:%M")
ST_UTC0=ST_UTC8-timedelta(hours=8)
ET_UTC0=ET_UTC8-timedelta(hours=8)
FileSet,times=DasFileRead(ST_UTC0,ET_UTC0,DataPath)
ST=times[0]
ET=times[-1]
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

#%%滤波


def bandpass_f(some_data,sample_rate,cutoff_f_low,cutoff_f_high,cutoff_order):
    #带通滤波
    data=[]
    cutoff_f0=cutoff_f_low/sample_rate*2 #截止频率10hz
    cutoff_f1=cutoff_f_high/sample_rate*2 
    #cutoff_order=cutoff_order   #配置滤波器的阶数
    if cutoff_f0==0:
        b,a = signal.butter(cutoff_order, cutoff_f1, 'low') 
    else:
        b,a = signal.butter(cutoff_order, [cutoff_f0,cutoff_f1], 'bandpass') 
    #分别对所有的通道滤波，滤波后的数据保存在some_data里
    for channel in range(0,some_data.shape[1]):
        data.append(signal.filtfilt(b,a, some_data[:,channel]))
        if (some_data.shape[1]-channel)%100==0:
            print('The rest channels for filtering:{}'.format(some_data.shape[1]-channel))
    return np.array(data).transpose()

FILTER=0
if FILTER==1:
    FILTER_Data = bandpass_f(Data, fs, 0,5,4) 
else:
    FILTER_Data=Data

#%%
#下采样（台站的采样率一般是100HZ以下）-----------------------------------------------------------------------------------------------------------------------------

DataCoordX, DataCoordy = FILTER_Data.shape
print('DataCoordX',DataCoordX,'DataCoordy',DataCoordy)
# 数据切片
t_start = 0
t_end = -1
t_interval = 1
TWindow = slice(t_start, t_end, t_interval)
c_start = 0
c_end = -1
c_interval = 1
CWindow = slice(c_start, c_end, c_interval)

DataSlice = FILTER_Data[TWindow, CWindow]

# 对some_data 进行采样，因为原始数据每秒采样1000，可以降为200个点以方便人为的判断
DownSampleRate = 100 #输入的采样数据为1秒1000个点，这里设置每秒采样的点数
FreqDownSample = DownSampleRate  #下采样频率
TDownSample = slice(0, DataCoordX, int(1000 / DownSampleRate))
DataDownSample = DataSlice[TDownSample, :]
print('DataDownSample',np.shape(DataDownSample))

#降采样后再滤波，节省计算量
fs=DownSampleRate
#FILTER_Data = bandpass_f(DataDownSample, fs,10,50,4) 
#高斯滤波
print('Filter!')
FILTER_Data=DataDownSample

#读取时间区间的船只数据
PosFile='pos_zf_gsd_1657814400_1660579199_742.csv'
StaticFile='static_zf_gsd_1657814400_1660579199_742.csv'
FiberBoatMessage=AISData(PosFile,StaticFile,ST,ET)

#%%画出振动时空图

Tstart = 1  #minute
Tend = -1   #minute
TimeWin = slice(Tstart*60*DownSampleRate, Tend*60*DownSampleRate, 1)
Cstart = 0
Cend= -1
Cwin = slice(Cstart,Cend,1)

ShowData=(FILTER_Data[TimeWin,Cwin])

CT=list(FiberBoatMessage['CrossTime'])
deltaCT=[]
for i in range(0,len(CT)):
    deltaCT.append(round((CT[i]-ST)/(ET-ST)*ShowData.shape[0]))

vline_indx = deltaCT

#利用AIS数据在DAS数据上标定船只通过的时间以及位置
Dist=list(FiberBoatMessage['disFromEnd/Km'])
#区域距离标定(3.33和120都是经验参数)
#RegionDistUpper=list(n_channels-(Dist+3.33)*1000/channel_spacing+120)
RegionDistUpper=[round(n_channels-(i+3.33)*1000/channel_spacing+200) for i in Dist]
RegionDistdown=[round(n_channels-(i+3.33)*1000/channel_spacing-200) for i in Dist]


deltaCTUpper=[i+round(ShowData.shape[0]*0.05) for i in vline_indx]
deltaCTdown=[i-round(ShowData.shape[0]*0.05) for i in vline_indx]


CT=list(FiberBoatMessage['Time_0'])
deltaCT=[]
for i in range(0,len(CT)):
    deltaCT.append(round((CT[i]-ST)/(ET-ST)*ShowData.shape[0]))


deltaCTdown=deltaCT

CT=list(FiberBoatMessage['Time_1'])
deltaCT=[]
for i in range(0,len(CT)):
    deltaCT.append(round((CT[i]-ST)/(ET-ST)*ShowData.shape[0]))
deltaCTUpper=deltaCT

fig1=plt.figure(dpi=1000)    
ax1 = fig1.add_subplot(1,1,1)
plt.imshow(((np.transpose(ShowData))), cmap="seismic", aspect='auto',origin='lower',vmin= 100,vmax=600) # cmap='seismic',, 

#计算坐标轴的刻度大小,合计10个时间戳(计算有问题，需要考虑数据的实际距离以及截断)
TimeTicks=10
TI=(ET-ST)/TimeTicks

xlabel=np.linspace(0,ShowData.shape[0],TimeTicks+1)
plt.xticks(xlabel,pd.date_range(ST.strftime("%Y%m%d %H%M"),ET.strftime("%Y%m%d %H%M"),freq=TI),rotation = 60)
ylabel=np.arange(0,n_channels,500)
plt.yticks(ylabel,np.round(ylabel*channel_spacing/1000))


#区域标定
for i in range(0,len(deltaCTdown)):
    plt.plot([deltaCTdown[i],deltaCTUpper[i],deltaCTUpper[i],deltaCTdown[i],deltaCTdown[i]],[RegionDistdown[i],RegionDistdown[i],RegionDistUpper[i],RegionDistUpper[i],RegionDistdown[i]],linewidth=1)

#plt.vlines(vline_indx, 0, 4000, colors='g', linestyles='dashed', label='垂直线')
plt.savefig("Space-time diagram.png")
plt.show()
print("Space-time diagram.png")
#%%
Ccord=ShowData.shape[0]
Tcord=ShowData.shape[1]
plt.figure(dpi=800)
plt.plot(np.arange(0,4094,1),ShowData[int(Ccord/2),:])
plt.xlabel('Channel')
plt.ylabel('Amplitude')
plt.savefig("Slice diagram.png")

# %%
