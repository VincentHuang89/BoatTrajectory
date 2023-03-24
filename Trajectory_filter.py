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
from RadonWake import PlotRadon,SpeedOnRadon,SNROnRadon,ValidationSpeedOnRadon,EnhanceResolution,CalculateKEnv
from Paperfigure import PlotK_KenvLine,PlotRadonInPaper,PlotSimulInDAS
from scipy import stats
#%%Params ---------------------------------------
#Save Params or not
SaveParams=0
# filter or not
FILTER=1
DownSampleRate = 50   #输入的采样数据为1秒1000个点，这里设置每秒采样的点数
# Showdata params
MINTIME=0.5   #0.5  0
MAXTIME=2   #2    3
MINCHANNEL=8.5   #8.5   7.8  
MAXCHANNEL=9.7 #Km  #9.7  10.5
#波线方向（左下到右上：0 （图像域的上半部分，deg：0-90），左上到右下：1，deg：90-180，展示所有：2）
WAVEDIRECT=1

# Z-score and threshold filter
threshold=1.5
#radon transfromation params
Tstart =2
Tend =2.5
Cstart =8.3
Cend= 9
CSTART=Cstart
CEND=Cend
Cstart=max(Cstart-MINCHANNEL,0)
Cend=min(Cend-MINCHANNEL,MAXCHANNEL)
#To anchor the ship or not in DAS figure with AIS data
PLOTANCHOR=0

#To anchor the dataslice based on Tstart,Tend,Cstart,Cend
PLOTREGION=1

#Denoise Radon-domain data or not 
DENOISE_RADON=1

#Enhance the resolution of DataSlice in channel and time dimension or not
SPACE_EN=0
TIME_EN=0
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
SHIP=0
MMSI=['413260090','413208430','413471740']
V=[13.86,13.82,12.83]
Angle=[119.86,58.22,123.26]
h=7

StartTime=["24/07/22 09:54","24/07/22 10:02","24/07/22 21:09"]
EndTime=["24/07/22 10:03","24/07/22 10:03","24/07/22 21:11"]

ST_UTC8=datetime.datetime.strptime(StartTime[SHIP], "%d/%m/%y %H:%M")
ET_UTC8=datetime.datetime.strptime(EndTime[SHIP], "%d/%m/%y %H:%M")
ST_UTC0=ST_UTC8-timedelta(hours=8)
ET_UTC0=ET_UTC8-timedelta(hours=8)
FileSet,times=DasFileRead(ST_UTC0,ET_UTC0,DataPath)
ST=times[0]
ET=times[-1]+timedelta(minutes=1)
ST=ST+timedelta(hours=8)
ET=ET+timedelta(hours=8)
print('Minutes of the DAS data is ',len(FileSet))

ST1=ST+timedelta(minutes=MINTIME)
if MAXTIME==-1:
    ET1=ET
else:
    ET1=ST+timedelta(minutes=MAXTIME)
print(ST1,ET1)


STW=ST+timedelta(minutes=Tstart)
if Tend==-1:
    ETW=ET
else:
    ETW=ST+timedelta(minutes=Tend)
print(STW,ETW)


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



#%%
'''
PlotDAS(ShowData,ST1,ET1,FiberBoatMessage,MINCHANNEL,MAXCHANNEL,RegionSliceX,RegionSliceY,channel_spacing,n_channels,PLOTANCHOR,PLOTREGION)  
#DASNR(ShowData)


#For Radon transformation, Slice partial data from ShowData-----------------------------
TimeWin = slice(int(Tstart*60*DownSampleRate), max(-1,int(Tend*60*DownSampleRate)), 1)
Cwin = slice(int(Cstart*1000/channel_spacing),max(-1,int(Cend*1000/channel_spacing)),1)
ShowDataSlice=(ShowData[TimeWin,Cwin])

#readjust the time dimension of the ShowDataSlice, to avoid the time-consuming calcualtion in the radon transformation caused by the too long time dimension
ShowDataSlice,ReDownSampleRate,channel_spacing_scaling=EnhanceResolution(ShowDataSlice,DownSampleRate,SPACE_EN,TIME_EN)
print('ReDownSampleRate',ReDownSampleRate)
ShowDataSlice=np.transpose((ShowDataSlice))

# Radon transformation and analysis
sinogram=PlotRadon(ShowDataSlice,DENOISE_RADON)
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
ShowDataSlice=ValidationSpeedOnRadon(speed,FILTER_Data,ReDownSampleRate,channel_spacing,MINTIME,MAXTIME,MINCHANNEL,MAXCHANNEL,Tstart,Tend,Cstart,Cend,threshold,WAVEDIRECT)

K_env,bias=CalculateKEnv(ShowDataSlice,channel_spacing,ReDownSampleRate)
Env_speed=K_env*channel_spacing*DownSampleRate

print('包络线速度',Env_speed,K_env)
#To validate the K-line and K-envelope-line in the same figure

#figures in the paper
PlotK_KenvLine(ShowDataSlice,speed,DownSampleRate,channel_spacing,WAVEDIRECT,K_env,bias,STW,ETW,CSTART)
PlotRadonInPaper(ShowDataSlice,channel_spacing,STW,ETW,CSTART)
'''

WLen_Scale=25   #20
Wbias=30  #63
Tbias=0.32  #0.3
h=7.246
angle=118.8
#angle=angle-90  #重新映射到光纤的角度
v=13.69
UpperBound=50
LowerBound=-10000
PlotSimulInDAS(DownSampleRate,v,h,angle,25,ShowData,ST1,ET1,MINCHANNEL,channel_spacing,WLen_Scale,Wbias,Tbias,UpperBound,LowerBound)
#Save setting Params to excel
if SaveParams==1:
    df=pd.DataFrame(np.array([FILTER,DownSampleRate,MINTIME,MAXTIME,MINCHANNEL,WAVEDIRECT,threshold,Tstart,Tend,CSTART,CEND,ST_UTC8,ET_UTC8,speed]).reshape(1,-1),columns=["FILTER","DownSampleRate","MINTIME","MAXTIME","MINCHANNEL","WAVEDIRECT","threshold","Tstart","Tend","Cstart","Cend","ST_UTC8","ET_UTC8","Speed of wave along the fiber"])

    if not os.path.exists('Analyze_DAS_Data.xlsx'):
        df.to_excel("Analyze_DAS_Data.xlsx")

    else: 
        res = pd.read_excel('Analyze_DAS_Data.xlsx',index_col=0)
        df = pd.concat([res,df],axis=0,ignore_index=1)
        df.to_excel("Analyze_DAS_Data.xlsx")
