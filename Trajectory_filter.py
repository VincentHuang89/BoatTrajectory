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
from AISData import AISData,AnchorShip
from DASFilter import bandpass_f,WeightMatrix
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
ST_UTC8=datetime.datetime.strptime("24/07/22 21:05", "%d/%m/%y %H:%M")
ET_UTC8=datetime.datetime.strptime("24/07/22 21:16", "%d/%m/%y %H:%M")
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
MINCHANNEL=0
MAXCHANNEL=-1
FILTER_Data=DataDownSample[:,slice(MINCHANNEL,MAXCHANNEL,1)]

#%%

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

MinValue=200
#ShowData=np.maximum(ShowData,MinValue)

#归一化
from sklearn import preprocessing
ShowData = preprocessing.scale(ShowData)


fig1=plt.figure(dpi=400,figsize=(13,10))    
ax1 = fig1.add_subplot(1,1,1)
#plt.imshow((np.abs(np.transpose(ShowData))), cmap="seismic", aspect='auto',origin='lower',vmin= 0,vmax=600) # cmap='seismic',, 
#plt.imshow(((np.transpose(ShowData))), cmap="bwr", aspect='auto',origin='lower',vmin= MinValue,vmax=700) # cmap='seismic',, 
plt.imshow(((np.transpose(ShowData))), cmap="bwr", aspect='auto',origin='lower',vmin=-3,vmax=3) # cmap='seismic',, 

#计算坐标轴的刻度大小,合计10个时间戳(计算有问题，需要考虑数据的实际距离以及截断)
TimeTicks=10
TI=(ET-ST)/TimeTicks
xlabel=np.linspace(0,ShowData.shape[0],TimeTicks+1)
plt.xticks(xlabel,pd.date_range(ST.strftime("%Y%m%d %H%M"),ET.strftime("%Y%m%d %H%M"),freq=TI),rotation = 60)
ylabel=np.arange(0,ShowData.shape[1],500)
plt.yticks(ylabel,np.round(ylabel*channel_spacing/1000))


#船只过光纤轨迹标定
PLOTANCHOR=0
if PLOTANCHOR==1:
    deltaCTUpper,deltaCTdown,RegionDistUpper,RegionDistdown=AnchorShip(FiberBoatMessage,MINCHANNEL,n_channels,channel_spacing,ST,ET,ShowData.shape[0])
    for i in range(0,len(deltaCTdown)):
        plt.plot([deltaCTdown[i],deltaCTUpper[i],deltaCTUpper[i],deltaCTdown[i],deltaCTdown[i]],[RegionDistdown[i],RegionDistdown[i],RegionDistUpper[i],RegionDistUpper[i],RegionDistdown[i]],linewidth=1)

#plt.vlines(vline_indx, 0, 4000, colors='g', linestyles='dashed', label='垂直线')
plt.savefig("Space-time diagram.png",bbox_inches = 'tight')
print("Space-time diagram.png")



# %%
'''

#%%1）速度滤波
L=2  #Km
t=4   #Min  #25.5,26.5,11.5,25.6
KernelSize=600
Speed=L*1000/(t*60)
K=Speed/(channel_spacing*DownSampleRate)
weight_matrix,Cdim,Tdim=WeightMatrix(K,KernelSize)
plt.figure()
plt.imshow(weight_matrix)
plt.plot(Tdim,Cdim)
#%%

Padding=int(KernelSize/2)
MaxTime=FILTER_Data.shape[0]
DataFilteredInSpeed=[]
for t0 in range(0,MaxTime-KernelSize):
    Tstart = t0  
    Tend = t0+KernelSize   
    TimeWin = slice(Tstart, Tend, 1)
    Cstart = 0
    Cend= -1
    Cwin = slice(Cstart,Cend,1)
    Data=((FILTER_Data[TimeWin,Cwin]))
    MaxChannel=Data.shape[1]


    DataFilteredInChannel=[]
    for c in range(0,MaxChannel-KernelSize):
        #利用channel维度，截取Kernel所对应的Ifmap数据
        Cstart = c
        Cend= c+KernelSize
        Cwin = slice(Cstart,Cend,1)
        DataKernel=Data[:,Cwin]
        #plt.figure()
        #plt.imshow(np.transpose(DataKernel), cmap="seismic",   aspect='auto',    origin='lower',vmin= 0,vmax=600)
        temp=0
        for t in Tdim:
            temp=temp+DataKernel[t][Cdim[t]]
        DataFilteredInChannel.append(temp/len(Tdim))
    DataFilteredInSpeed.append(DataFilteredInChannel)
    if t0%100==0:
        print('Filter has been executed '+str(t0)+' times/'+str(MaxTime))
DataFilteredInSpeed=np.array(DataFilteredInSpeed)

# %%
#plt.figure()
#plt.imshow(weight_matrix)
plt.figure(dpi=800)
plt.imshow((np.transpose(DataFilteredInSpeed)), cmap="seismic",   aspect='auto',    origin='lower',vmin= -200,vmax=200)
plt.savefig('SpeedFiltered.png')
# %%
'''