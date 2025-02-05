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
import matplotlib
import scipy.signal as signal
from scipy.fftpack import fft,ifft
from sklearn.decomposition import FastICA
import datetime
from datetime import timedelta
from DASFileRead import DasFileRead
from tdms_reader_1 import *
from AISData import AISData,AnchorShip,Index_AIS_csv_by_DAS_Data
from DASFilter import bandpass_f,DWT_all_channels,DataDiff,WeightMatrix,PlotDAS,DASNR,DASInterpolate,Dup_Spat_Dim,Dup_Time_Dim,Pos_AnchorRegion,Cal_pos_s,Cal_pos_t,Cal_region_all_pos,Params_ShowData_3slice
import skimage 
from skimage.transform import radon
from RadonWake import PlotRadon,SpeedOnRadon,SNROnRadon,ValidationSpeedOnRadon,EnhanceResolution,CalculateKEnv
from Paperfigure import PlotK_KenvLine,PlotRadonInPaper,PlotSimulInDAS
from scipy import stats
plt.rcParams['font.sans-serif']=['SimHei'] #用来正常显示中文标签
plt.rcParams['axes.unicode_minus']=False #用来正常显示负号

os.chdir('/home/huangwj/DAS/BoatTrajectory')
print('当前工作路径是：' + os.getcwd())  # 显示当前路径

#%%Params ---------------------------------------
SHIP=0
MMSI=['413260090','413208430','413471740','413226010','413231470','413260090']
DataPath='/home/huangwj/DAS/BoatTrajectory/DataforAIS'
StartTime=["24/07/22 09:54","24/07/22 10:01","24/07/22 21:09",'10/09/21 13:17','06/09/21 11:56','09/09/21 9:33']
EndTime=["24/07/22 10:00","24/07/22 10:04","24/07/22 21:11",'10/09/21 13:22','06/09/21 11:58','09/09/21 9:37']
V=[13.86,13.82,12.83,15.07,16.15,14.33]
Angle=[119.86,58.22,123.26,88.39,90.01,46.63]
h=7

#定义船行波尖峰时刻和位置
Ship_relocate = [('24/07/22 09:54:40',9.3),
                 ('24/07/22 10:02:10',11.25),
                 ('24/07/22 21:09:40',11.1),
                 ('10/09/21 13:18:30',12.55),
                 ('06/09/21 11:56:10',11.6),
                 ('09/09/21 9:34:00',9.58)]

#定义数据切片左上顶点到船行波尖峰的时间差距以及位置差距；定义切片时间长度以及距离长度
#Bias_region_shipwave_top = [(5/6,-0.15),(5/6,-0.15),(5/6,-0.15),(5/6,-0.15),(0.5,0.5),(0.5,-0.10)]  #bias_t(min),bias_s(Km)
#len_region = [(0.5,0.6),(0.5,0.6),(0.5,0.6),(0.5,0.6),(0.5,0.5),(0.3,0.2)]  #min,Km

Bias_region_shipwave_top_set = {'k2':[(1.33,1.2),(0.3,0.45),(3/6,0.5),(2/6,0.3),(0.5,0.5),(0.2,0.3 )],
                                'k1':[(1.33,-0.25),(0.3,-0.1),(3/6,-0.08),(2/6,-0.08),(0.5,-0.1),(0.2,-0.15)],
                                'ke':[(0.5,-0.1),(0.3,0.45),(3/6,-0.08),(2/6,-0.07),(0.5,-0.1),(0.2,0.3)] }
len_region_set = {'k2':[(0.5,0.6),(0.5,0.4),(0.5,0.4),(0.5,0.2),(0.3,0.4),(0.3,0.2)],
                  'k1':[(0.5,0.6),(0.5,0.4),(0.5,0.4),(0.5,0.2),(0.3,0.4),(0.3,0.2)],
                  'ke':[(0.75,0.8),(0.75,0.4),(0.75,0.6),(0.8,0.6),(0.5,0.5),(1.6,0.2)]}  #min,Km

Lateral_Wave =['K1','K2','K1','K1','K1']  #用于分析不同单边散波
Lateral_Wave_selected= Lateral_Wave[SHIP]
#Save Params or not
SaveParams=0
# filter or not
FILTER = 1
SamplesPerSec = 10 #每秒采样点

# Showdata params
MINTIME = [0,0,0,0,0,0]   #0.5  0
MAXTIME = [-1,-1,-1,-1,2.5,-1]   #2    3
MINCHANNEL=[8.5,10,10,11,10,9]   #8.5   7.8  
MAXCHANNEL=[10.5,12.7,12.7,12.9,12.5,10.5] #Km  #9.7  10.5
MINTIME=MINTIME[SHIP]
MAXTIME=MAXTIME[SHIP]
MINCHANNEL=MINCHANNEL[SHIP]
MAXCHANNEL=MAXCHANNEL[SHIP]



# Z-score and threshold filter
threshold = [1.2,1,0.5,1.5,1.2,0] #1.6
#threshold = [0,0,0.2,0.2,0.5,0] #1.6
threshold=threshold[SHIP]

#radon transfromation params   (!!!先重定位船行波尖峰具体的位置，然后在定位下述区域，避免调整数据，导致区域偏移)
Tstart = [1,1,1,1,3,1]
Tend = [1.5,1.5,1.5,1.5,3.5,3]
Cstart = [12,12,12,12,10.85,6]
Cend= [1.5,1.5,1.5,1.5,11.45,12]

Tstart=Tstart[SHIP]
Tend=Tend[SHIP]
Cstart=Cstart[SHIP]
Cend=Cend[SHIP]

CSTART=Cstart
CEND=Cend
Cstart=max(Cstart-MINCHANNEL,0)
Cend=min(Cend-MINCHANNEL,MAXCHANNEL)
#To anchor the ship or not in DAS figure with AIS data
PLOTANCHOR = 0

#To anchor the dataslice based on Tstart,Tend,Cstart,Cend
PLOTREGION = 1

#Denoise Radon-domain data or not 
DENOISE_RADON = 1

#Enhance the resolution of DataSlice in channel and time dimension or not
SPACE_EN = 0
TIME_EN = 0
#---------------------------------------------------

mpl.rcParams['font.sans-serif'] = ['SimHei']  # 显示中文
mpl.rcParams['axes.unicode_minus'] = False  # 显示负号



#Read AIS data .csv files
AIS_file_list = Index_AIS_csv_by_DAS_Data(datetime.datetime.strptime(StartTime[SHIP], "%d/%m/%y %H:%M"))
print(AIS_file_list)
FiberBoatMessage = pd.read_csv(AIS_file_list[0])
FiberBoatMessage['CrossTime'] = pd.to_datetime(FiberBoatMessage['CrossTime'],format="%Y-%m-%d %H:%M:%S")
FiberBoatMessage['Time_1'] = pd.to_datetime(FiberBoatMessage['Time_1'],format="%Y-%m-%d %H:%M:%S")
FiberBoatMessage['Time_0'] = pd.to_datetime(FiberBoatMessage['Time_0'],format="%Y-%m-%d %H:%M:%S")


#锚定船行波出现的时间和位置

ST_UTC8=datetime.datetime.strptime(StartTime[SHIP], "%d/%m/%y %H:%M")
ET_UTC8=datetime.datetime.strptime(EndTime[SHIP], "%d/%m/%y %H:%M")
ST_UTC0=ST_UTC8-timedelta(hours=8)
ET_UTC0=ET_UTC8-timedelta(hours=8)
print(DataPath)
FileSet,times=DasFileRead(ST_UTC0,ET_UTC0,DataPath)
print(FileSet)
print(times)
ST=times[0]
ET=times[-1]+timedelta(minutes=1)
ST=ST+timedelta(hours=8)
ET=ET+timedelta(hours=8)
print('Minutes of the DAS data is ',len(FileSet))
print(FileSet)

ST1=ST+timedelta(minutes=MINTIME)  #在ST-ET区间中按照MINTIME和MAXTIME重新选择数据
if MAXTIME==-1:
    ET1=ET
else:
    ET1=ST+timedelta(minutes=MAXTIME)

print(ST1,ET1)



#%%Read DAS data files
start = time.time()

all_data_downsampled = pd.DataFrame()
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
        print(f'Zero_offset is:{zero_offset}')
        print('Channel_spacing is:{:.4f}m'.format(channel_spacing))
        print('Number of channels in file: {0}'.format(n_channels))
        print('Time samples in file: {0}'.format(tdms.channel_length))
        print('Sampling frequency (Hz): {0}'.format(fs))
        print(f'Total length: {np.max(depth)}m')
        down_sample_factor = int(props.get('SamplingFrequency[Hz]')/SamplesPerSec)
        first_channel = 0  
        last_channel = n_channels - 1  
        first_time_sample = 00 
        last_time_sample = tdms.channel_length - 1 
        #some_data = tdms.get_data(first_channel, last_channel, first_time_sample, last_time_sample)   
        some_data = pd.DataFrame(tdms.get_data(first_channel, last_channel, first_time_sample, last_time_sample))
        data_downsampled = some_data.groupby(some_data.index // down_sample_factor).mean()  # 平均采样
        del some_data  # 释放内存
        all_data_downsampled = pd.concat([all_data_downsampled, data_downsampled], axis=0)
    else:
        tdms = TdmsReader(FilePath)
        some_data = pd.DataFrame(tdms.get_data(first_channel, last_channel, first_time_sample, last_time_sample))        
        data_downsampled = some_data.groupby(some_data.index // down_sample_factor).mean()  # 平均采样
        del some_data  # 释放内存
        all_data_downsampled = pd.concat([all_data_downsampled, data_downsampled], axis=0)
    print(FilePath)

#%% Mean Sampling
#Mean sampling

print('Size of downsampled data: {0}'.format(all_data_downsampled.shape))
DataDownSample = np.array(all_data_downsampled)

REGION = Params_ShowData_3slice(Ship_relocate,ST1,MINCHANNEL,Bias_region_shipwave_top_set,len_region_set,SamplesPerSec,channel_spacing,SHIP)
print(REGION)

STW=ST+timedelta(minutes=Tstart)
if Tend==-1:
    ETW=ET
else:
    ETW=ST+timedelta(minutes=Tend)



#Get the differential data in the time dimension
#FILTER_Data=DataDiff(FILTER_Data,1)

#%%

#Obtain the data according to the Time and Channel Parameters from the downsampled data, to show the DAS image
TimeWin = slice(int(MINTIME*60*SamplesPerSec), max(-1,int(MAXTIME*60*SamplesPerSec)), 1)
ShowData = DataDownSample[TimeWin,slice(int(MINCHANNEL*1000/channel_spacing),int(MAXCHANNEL*1000/channel_spacing),1)]

#%%

ShowData = stats.zscore(ShowData,axis = 0)
ShowData_raw = ShowData

#%%Filter the raw data
if FILTER==1:
    #FILTER_Data = bandpass_f(DataDownSample, fs, 0,5,4) 
    ShowData = DWT_all_channels(ShowData,wavelet = 'db9',level = 7)
    print('Filter!')

ShowData = stats.zscore(ShowData, axis=None)
STD=threshold*np.std(ShowData,ddof=1)
ShowData=((ShowData>STD)|(ShowData<-STD))*ShowData


PlotDAS(ShowData,ST1,ET1,FiberBoatMessage,MINCHANNEL,MAXCHANNEL,REGION,channel_spacing,n_channels,zero_offset,PLOTANCHOR,PLOTREGION,SHIP)  



print(ShowData.shape)

#For Radon transformation, Slice partial data from ShowData-----------------------------
Lateral_Wave = 'K1'
Region_all_pos = REGION[Lateral_Wave]
TimeWin = slice(Region_all_pos[0][0],Region_all_pos[1][0],1)
Cwin = slice(Region_all_pos[3][1],Region_all_pos[0][1],1)
ShowDataSlice=(ShowData[TimeWin,Cwin])

#readjust the time dimension of the ShowDataSlice, to avoid the time-consuming calcualtion in the radon transformation caused by the too long time dimension
ShowDataSlice,ReDownSampleRate,channel_spacing_scaling=EnhanceResolution(ShowDataSlice,SamplesPerSec,SPACE_EN,TIME_EN)
print('ReDownSampleRate',ReDownSampleRate)
#ReDownSampleRate=DownSampleRate
#channel_spacing_scaling=1
ShowDataSlice=np.transpose((ShowDataSlice))

# Radon transformation and analysis
sinogram=PlotRadon(ShowDataSlice,DENOISE_RADON)
print("Radon done!")
deg_mean,deg_median,s_mean,s_median=SpeedOnRadon(sinogram,max(ShowDataSlice.shape),channel_spacing,ReDownSampleRate,channel_spacing_scaling,Lateral_Wave)

res=[]
res.append((ReDownSampleRate,deg_mean,deg_median,s_mean,s_median))
RES=pd.DataFrame(res,columns=['DownSampleRate',"deg_mean","deg_median","s_mean",'s_median'])
#RES.to_excel('res.xlsx')
RES.dropna(inplace=True)
s_mean = RES['s_mean']
speed_1=np.mean(s_mean)
print('Estimated ship speed(K1): ',speed_1)


#To validate the accuracy of the estimated speed, plot the line according to the estimated speed in the ShowDataSlice image.
#ShowDataSlice=ValidationSpeedOnRadon(speed,FILTER_Data,ReDownSampleRate,channel_spacing,MINTIME,MAXTIME,MINCHANNEL,MAXCHANNEL,Region_all_pos,threshold,WAVEDIRECT)



#For Radon transformation, Slice partial data from ShowData-----------------------------
Lateral_Wave = 'K2'
Region_all_pos = REGION[Lateral_Wave]
TimeWin = slice(Region_all_pos[0][0],Region_all_pos[1][0],1)
Cwin = slice(Region_all_pos[3][1],Region_all_pos[0][1],1)
ShowDataSlice=(ShowData[TimeWin,Cwin])

#readjust the time dimension of the ShowDataSlice, to avoid the time-consuming calcualtion in the radon transformation caused by the too long time dimension
ShowDataSlice,ReDownSampleRate,channel_spacing_scaling=EnhanceResolution(ShowDataSlice,SamplesPerSec,SPACE_EN,TIME_EN)
print('ReDownSampleRate',ReDownSampleRate)
#ReDownSampleRate=DownSampleRate
#channel_spacing_scaling=1
ShowDataSlice=np.transpose((ShowDataSlice))

# Radon transformation and analysis
sinogram=PlotRadon(ShowDataSlice,DENOISE_RADON)
print("Radon done!")
deg_mean,deg_median,s_mean,s_median=SpeedOnRadon(sinogram,max(ShowDataSlice.shape),channel_spacing,ReDownSampleRate,channel_spacing_scaling,Lateral_Wave)

res=[]
res.append((ReDownSampleRate,deg_mean,deg_median,s_mean,s_median))
RES=pd.DataFrame(res,columns=['DownSampleRate',"deg_mean","deg_median","s_mean",'s_median'])
#RES.to_excel('res.xlsx')
RES.dropna(inplace=True)
s_mean = RES['s_mean']
speed_2=np.mean(s_mean)
print('Estimated ship speed(K2): ',speed_2)

lateral_speed = {'K1':speed_1,'K2':speed_2}

#包络线-------------

Region_all_pos = REGION['Ke']
TimeWin = slice(Region_all_pos[0][0],Region_all_pos[1][0],1)
Cwin = slice(Region_all_pos[3][1],Region_all_pos[0][1],1)
ShowDataSlice=np.transpose(ShowData[TimeWin,Cwin])

K_env,bias=CalculateKEnv(ShowDataSlice,channel_spacing,ReDownSampleRate,Lateral_Wave_selected)
Env_speed=K_env*channel_spacing*SamplesPerSec
print('包络线速度',Env_speed,K_env)
#To validate the K-line and K-envelope-line in the same figure

#figures in the paper
PlotK_KenvLine(ShowDataSlice,lateral_speed[Lateral_Wave_selected],SamplesPerSec,channel_spacing,Lateral_Wave_selected,K_env,bias,STW,ETW,CSTART,SHIP)
PlotRadonInPaper(ShowDataSlice,channel_spacing,STW,ETW,CSTART,SHIP)



'''
OccurTime=["24/07/22 09:54:35"]
OccurDist=[9.25]
OCCURTIME=datetime.datetime.strptime(OccurTime[SHIP], "%d/%m/%y %H:%M:%S")
OCCURDIST=OccurDist[SHIP]
WLen_Scale=15   #20
Wbias=OCCURDIST-MINCHANNEL
Tbias=(OCCURTIME-ST1)/timedelta(seconds=1)
h=7.246
angle=118.8
#angle=angle-90  #重新映射到光纤的角度
v=13.69
UpperBound=50
LowerBound=-10000
PlotSimulInDAS(DownSampleRate,v,h,angle,25,ShowData,ST1,ET1,MINCHANNEL,channel_spacing,WLen_Scale,Wbias,Tbias,UpperBound,LowerBound)
'''


REGION['Wave_peak']

OccurTime,OccurDist=Ship_relocate[SHIP]

OCCURTIME=datetime.datetime.strptime(OccurTime, "%d/%m/%y %H:%M:%S")
OCCURDIST=OccurDist
WLen_Scale=10  #20
Wbias=OCCURDIST-MINCHANNEL
Tbias=(OCCURTIME-ST1)/timedelta(seconds=1)
h=9.991394444
angle=100.4360294
#angle=angle-90  #重新映射到光纤的角度
v=15.68082052
UpperBound=500
LowerBound=-500
PlotSimulInDAS(REGION,SamplesPerSec,v,h,angle,80,ShowData_raw,ST1,ET1,MINCHANNEL,channel_spacing,WLen_Scale,Wbias,Tbias,UpperBound,LowerBound,SHIP)


#Save setting Params to excel
if SaveParams==1:
    df=pd.DataFrame(np.array([FILTER,SamplesPerSec,MINTIME,MAXTIME,MINCHANNEL,threshold,Tstart,Tend,CSTART,CEND,ST_UTC8,ET_UTC8,(speed_1,speed_2)]).reshape(1,-1),columns=["FILTER","DownSampleRate","MINTIME","MAXTIME","MINCHANNEL","threshold","Tstart","Tend","Cstart","Cend","ST_UTC8","ET_UTC8","Speed of wave along the fiber"])

    if not os.path.exists('Analyze_DAS_Data.xlsx'):
        df.to_excel("Analyze_DAS_Data.xlsx")

    else: 
        res = pd.read_excel('Analyze_DAS_Data.xlsx',index_col=0)
        df = pd.concat([res,df],axis=0,ignore_index=1)
        df.to_excel("Analyze_DAS_Data.xlsx")
