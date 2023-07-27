#研究海浪噪声的分布
#%%
from cmath import exp
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

#mpl.use('tKAgg')
import matplotlib.pyplot as plt
import scipy.signal as signal
from scipy.fftpack import fft,ifft
from sklearn import preprocessing
from sklearn.decomposition import FastICA
from distfit import distfit 
from fitter import Fitter
from tdms_reader_1 import *


start = time.time()

mpl.rcParams['font.sans-serif'] = ['SimHei']  # 显示中文
mpl.rcParams['axes.unicode_minus'] = False  # 显示负号
# from func import *


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
select = 1
if select == 0 :   #三角岛220708 14:00 ~14:30 UTC+8
    fileRange = [100,117]  #文件的最大值是
    dataPath = '/home/huangwj/DAS/BoatTrajectory/DataforAIS/sanjiao-1k-4m_UTC_20220808_010038.863.tdms'
    dataPath1 = ['/home/huangwj/DAS/BoatTrajectory/DataforAIS/sanjiao-1k-4m_UTC_20220808_0','38.863.tdms']
elif select ==1:
    fileRange = [436,440]  #文件的最大值是 [436,449] 
    dataPath = '/home/huangwj/DAS/BoatTrajectory/DataforAIS/sanjiao-1k-4m_UTC_20220808_044038.863.tdms'
    dataPath1 = ['/home/huangwj/DAS/BoatTrajectory/DataforAIS/sanjiao-1k-4m_UTC_20220808_0','38.863.tdms']
elif select ==2:
    fileRange = [204,222]  #文件的最大值是
    dataPath = '/home/huangwj/DAS/BoatTrajectory/DataforAIS/sanjiao-1k-4m_UTC_20220808_120438.863.tdms'
    dataPath1 = ['/home/huangwj/DAS/BoatTrajectory/DataforAIS/sanjiao-1k-4m_UTC_20220808_1','38.863.tdms']
elif select ==3:
    fileRange = [621,660]  #文件的最大值是
    dataPath = '/home/huangwj/DAS/BoatTrajectory/DataforAIS/sanjiao-1k-4m_UTC_20220811_062138.863.tdms'
    dataPath1 = ['/home/huangwj/DAS/BoatTrajectory/DataforAIS/sanjiao-1k-4m_UTC_20220811_0','38.863.tdms']
else:
    fileRange = [529,560]  #文件的最大值是
    dataPath = '/home/huangwj/DAS/BoatTrajectory/DataforAIS/sanjiao-1k-4m_UTC_20220811_152938.863.tdms'
    dataPath1 = ['/home/huangwj/DAS/BoatTrajectory/DataforAIS/sanjiao-1k-4m_UTC_20220811_1','38.863.tdms']

# 设定工作路径
path = '/home/huangwj/DAS/BoatTrajectory/'
os.chdir(path)
print('当前工作路径是：' + os.getcwd())  # 显示当前路径
# 设定提取数据文件
# file_path = 'test4-1k-2m-rcv255-25ns_UTC_20210730_035343.492.tdms'
file_path = dataPath
print('File: {0}'.format(file_path))

# 读取文件信号，获取zero offset,SpatialResolution,n_channels,SamplingFrequency
tdms = TdmsReader(file_path)
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

# 截取文件中部分数据，时间开始和结束，位置开始和结束。
first_channel = 0  # 选择开始的位置
# first_channel =int((1700-zero_offset)/channel_spacing)
last_channel = n_channels - 1  # 选择结束位置,光纤距离= channel*SpatialResolution（m)
# last_channel= int((1800-zero_offset)/channel_spacing)

first_time_sample = 00  # 选择开始时间，单位ms
last_time_sample = tdms.channel_length - 1  # 30000-1 #选择截止时间，单位ms
# last_time_sample =30000

some_data = tdms.get_data(first_channel, last_channel, first_time_sample, last_time_sample)
Data = some_data


# 数据汇总--------------------------------------------------------------------------------------------------------------------------------
for i in range(fileRange[0]+1,fileRange[1],1):  #文件的最大值是1913
    file_path = dataPath1[0]+str(i)+dataPath1[1]
    print(file_path)
    tdms = TdmsReader(file_path)
    some_data1 = tdms.get_data(first_channel, last_channel, first_time_sample, last_time_sample)
    Data = np.concatenate((Data, some_data1), axis=0)


#下采样（台站的采样率一般是100HZ以下）-----------------------------------------------------------------------------------------------------------------------------
#%%
DataCoordX, DataCoordy = Data.shape
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

DataSlice = Data[TWindow, CWindow]

# 对some_data 进行采样，因为原始数据每秒采样1000，可以降为200个点以方便人为的判断
DownSampleRate = 10 #输入的采样数据为1秒1000个点，这里设置每秒采样的点数
FreqDownSample = DownSampleRate  #下采样频率
TDownSample = slice(0, DataCoordX, int(1000 / DownSampleRate))
DataDownSample = DataSlice[TDownSample, :]

#%%
Channel=2000
dist=distfit(todf=True)
dist.fit_transform(DataDownSample[:,Channel])
dist.plot()
#%%
f=Fitter(DataDownSample[:,Channel])
f.fit()
f.summary()
f.get_best()


# %%
Channles = slice(1000,2500,1)
Time = 1200
dist=distfit(distr=['norm','t','rayleigh'],todf=True)
dist.fit_transform(DataDownSample[Time,Channles])
dist.plot()

#%%
f=Fitter(DataDownSample[Time,Channles],distributions=['gamma', 'rayleigh', 'uniform','norm'], timeout =10000)
f.fit()
f.summary()
f.get_best()
#%%
f.fitted_param['norm']
f.hist()
# %%
f.plot_pdf(Nbest=5, lw=2, method='sumsquare_error')
# %%
