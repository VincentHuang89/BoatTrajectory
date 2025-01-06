#读取DAS数据保存在/home/huangwj/DAS/BoatTrajectory/Data
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


from tdms_reader_1 import *


start = time.time()

mpl.rcParams['font.sans-serif'] = ['SimHei']  # 显示中文
mpl.rcParams['axes.unicode_minus'] = False  # 显示负号
# from func import *


#选择数据路径
select = 1
if select == 1 :   #三角岛220708 14:00 ~14:30 UTC+8
    fileRange = [600,631]  #文件的最大值是630
    dataPath = '/home/huangwj/DAS/BoatTrajectory/Data/sanjiao-1k-4m_UTC_20220708_060042.090.tdms'
    dataPath1 = ['/home/huangwj/DAS/BoatTrajectory/Data/sanjiao-1k-4m_UTC_20220708_0','42.090.tdms']
else:
    fileRange = [38,48,1]  #文件的最大值是47
    dataPath = '/home/huangwj/DAS/DAS/Dark fiber/822-1839阳西1.9级地震/test3-1k-4m-rcv255-pls50ns_UTC_20210822_183836.917.tdms'
    dataPath1 = ['/home/huangwj/DAS/DAS/Dark fiber/822-1839阳西1.9级地震/test3-1k-4m-rcv255-pls50ns_UTC_20210822_18','36.917.tdms']

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
for i in range(fileRange[0],fileRange[1],1):  #文件的最大值是1913
    file_path = dataPath1[0]+str(i)+dataPath1[1]
    print(file_path)
    tdms = TdmsReader(file_path)
    some_data1 = tdms.get_data(first_channel, last_channel, first_time_sample, last_time_sample)
    Data = np.concatenate((Data, some_data1), axis=0)

# Data =np.diff(np.diff(Data,axis=0))

#%%滤波


def bandpass_f(some_data,sample_rate,cutoff_f_low,cutoff_f_high,cutoff_order):
    #带通滤波
    data=[]
    cutoff_f0=cutoff_f_low/sample_rate*2 #截止频率10hz
    cutoff_f1=cutoff_f_high/sample_rate*2 #截止频率10hz
    #cutoff_order=cutoff_order   #配置滤波器的阶数
    b,a = signal.butter(cutoff_order, [cutoff_f0,cutoff_f1], 'bandpass') 
    #分别对所有的通道滤波，滤波后的数据保存在some_data里
    for channel in range(0,some_data.shape[1]):
        data.append(signal.filtfilt(b,a, some_data[:,channel]))
        if (some_data.shape[1]-channel)%100==0:
            print('The rest channels for filtering:{}'.format(some_data.shape[1]-channel))
    return np.array(data).transpose()


#FILTER_Data = bandpass_f(Data, fs, 1,100,4)  #滤波器阶数不能设置太高


FILTER_Data=Data


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
DownSampleRate = 100  #输入的采样数据为1秒1000个点，这里设置每秒采样的点数
FreqDownSample = DownSampleRate  #下采样频率
TDownSample = slice(0, DataCoordX, int(1000 / DownSampleRate))
DataDownSample = DataSlice[TDownSample, :]
print('DataDownSample',np.shape(DataDownSample))


#%%
start_time = 0
end_time = 60
first_time_sample = int(start_time*fs)  #0#10000   #选择开始时间，单位ms
last_time_sample = int(end_time*fs)#最长的时间长度不可超过last_channel = n_channels-1

plot_derate=10   #绘制的降采率
    # plot_data=some_data[0::plot_derate,:]
extent=[first_time_sample/fs,Data[0::plot_derate,:].shape[0]/fs*plot_derate,depth[first_channel],depth[last_channel] ]
print(extent)
plt.figure()
vmin=-np.std(Data)*5
vmax=np.std(Data)*5
    #设置x y轴的值的范围
    #cmap为图片显示的颜色，autumn:红-橙-黄 gray:黑-白 hot:黑-红-黄-白 jet:蓝-青-黄-红(prefer)
print(vmin,vmax)
img1 = plt.imshow(Data.transpose(),cmap='seismic', aspect='auto',vmin=vmin,vmax=vmax,extent=extent,origin='lower')
plt.xlabel('Time (s)',fontsize=20)
plt.ylabel('Distance (m)',fontsize=20)
plt.show()
plt.savefig('Space-time.png')

print('Space-time.pnd Done')




'''


#%%画出振动时空图
Tstart = 0  #minute
Tend = 10   #minute
TimeWin = slice(Tstart*60*DownSampleRate, Tend*60*DownSampleRate, 1)
Cstart = 1700
Cend= 2700
Cwin = slice(Cstart,Cend,1)

# print(some_data)
plt.figure()
#plt.xlabel("Time"+"/interval "+str(1/DownSampleRate*1000)+" ms")
#plt.ylabel("Channel")
#plt.title(" Space-time diagram \n(Channel is from " + str(c_start) + " to " + str(c_end) + " (Space interval is " + str(c_interval) + "), \n Time is from " + str(t_start) + ":" + str(t_end) + " (Time interval is " + str(1/DownSampleRate*1000) + " ms))")
plt.imshow((np.abs(np.transpose(DataDownSample[TimeWin,Cwin]))), cmap='jet',vmin= 0,vmax=700) #, norm=mpl.colors.Normalize(0, 255)
#plt.pcolormesh(xx,yy,np.abs(np.transpose(DataDownSample)),shading= "auto",cmap = "jet")
plt.colorbar()
plt.savefig("Space-time diagram.png")
plt.show()
print("Space-time diagram.png")


#%%参看单通道的数据

Time_int = slice(0,60000,1)
Channel = 1200
plt.figure()
plt.plot(DataDownSample[Time_int,Channel])
plt.xlabel("time")
plt.ylabel("channel")
plt.savefig("test.png")


# %%
'''