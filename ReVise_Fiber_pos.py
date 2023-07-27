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
from DASFilter import bandpass_f,DataDiff,WeightMatrix,PlotDAS,DASNR,DASInterpolate,Dup_Spat_Dim,Dup_Time_Dim,Pos_AnchorRegion,Cal_pos_s,Cal_pos_t,Cal_region_all_pos
import skimage 
from skimage.transform import radon
from RadonWake import PlotRadon,SpeedOnRadon,SNROnRadon,ValidationSpeedOnRadon,EnhanceResolution,CalculateKEnv
from Paperfigure import PlotK_KenvLine,PlotRadonInPaper,PlotSimulInDAS
import download_from_FTP as FTP
import glob
import fnmatch
from scipy import stats

#读取筛选后的船只数据
filename = '/home/huangwj/DAS/BoatTrajectory/AIS_SHIP_DATA/SHIP_selected.csv'

Ship_filterd_pd = pd.read_csv(filename)
print(Ship_filterd_pd.columns)

Ship_filterd_pd['CrossTime'] = pd.to_datetime(Ship_filterd_pd['CrossTime'])


image_path = r'/home/huangwj/DAS/BoatTrajectory/Revise_Fiber_pos_Image'
DAS_temp_path = r'/home/huangwj/DAS/BoatTrajectory/DataforAIS'  # 存放DAS tdms文件的临时文件夹
ship_wave_fig_dir = r'/home/huangwj/DAS/DataSet/dataset/ship_raw'  # DAS图片数据集文件夹
FTP_dir = '/NASShare4/sanjiao_guishan0904_1117/'
ftp = FTP.create_FTP('172.18.218.194', 'hzd', 'optical503', FTP_dir)  # 创建ftp对象并登录选择指定文件夹
DAS_temp_files = glob.glob(os.path.join(DAS_temp_path, '*.tdms'))  # 获取DAS tdms临时文件列表

file_list = ftp.nlst()

DownSampleRate=100
MAX_DISTANCE = 12.9    # 岸上的最大距离(km)
MIN_DISTANCE = 2.7    # 岸上的最小距离(km)

TIME_left=-3
TIME_right = 3
for i in range(0,len(Ship_filterd_pd)):
    image_files = []
    image_filename = os.path.join(image_path,str(Ship_filterd_pd['MMSI'].iloc[i])+'_'+datetime.datetime.strftime(Ship_filterd_pd['CrossTime'].iloc[i],"%Y%m%d_%H%M")+'_'+str(round(Ship_filterd_pd['disFromEnd/Km'].iloc[i],2))+'.png')

    for j in range(TIME_left,TIME_right):
        datetime_str= datetime.datetime.strftime(Ship_filterd_pd['CrossTime'].iloc[i]-timedelta(hours=8)+timedelta(minutes=j),"%Y%m%d_%H%M")
        #print(Ship_filterd_pd['CrossTime'].iloc[i],datetime_str)
        try:
            for file in file_list:
                if fnmatch.fnmatch(file, f'*{datetime_str}*.tdms'):   # 通配符匹配文件名
                    print('Downloading:', file)
                    image_files.append(file)
                    local_file = os.path.join(DAS_temp_path, file)
                    FTP.download_file(ftp, file, local_file)
        except:
            print('下载被中断！')
    if len(image_files)<(TIME_right-TIME_left):
        print(f'{image_filename}的tdms 文件不完整，跳过')
        continue
    all_data_downsampled = pd.DataFrame()
    for k in range(0,len(image_files)):
        FilePath= os.path.join(DAS_temp_path,image_files[k])

        if k==0:
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
            down_sample_factor = fs/DownSampleRate
            first_channel = int(MIN_DISTANCE*1000/channel_spacing)  
            last_channel = int(MAX_DISTANCE*1000/channel_spacing)   
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

    #读取单条记录所需的数据，并画图
    ST = Ship_filterd_pd['CrossTime'].iloc[i]+timedelta(minutes=TIME_left)
    ET = Ship_filterd_pd['CrossTime'].iloc[i]+timedelta(minutes=TIME_right)
    

    ShowData= np.array(all_data_downsampled)
    ShowData = stats.zscore(ShowData, axis=None)

    Pass_time_dim = ShowData.shape[0]*abs(timedelta(minutes=TIME_left).total_seconds()/timedelta(minutes=(TIME_right-TIME_left)).total_seconds())

    fig1=plt.figure(dpi=400,figsize=(15,10))    
    #ax1 = fig1.add_subplot(1,1,1) 
    plt.imshow(np.transpose(ShowData), cmap="bwr", aspect='auto',origin='lower',vmin=-3,vmax=3,extent=[0, ShowData.shape[0], 0, ShowData.shape[1]]) # cmap=''bwr,
    plt.colorbar()
    TimeTicks = 7
    xlabel = np.linspace(0, ShowData.shape[0] - 1, TimeTicks)
    xtick_labels = pd.date_range(ST.strftime("%Y%m%d %H%M%S"), ET.strftime("%Y%m%d %H%M%S"), periods=TimeTicks).strftime('%H:%M:%S')
    TimeTicks=7
    xlabel=np.linspace(0,ShowData.shape[0]-1,TimeTicks)
    plt.xticks(xlabel,xtick_labels,rotation = 0,size=15)
    ylabel=np.linspace(0,ShowData.shape[1],5)
    plt.xlabel("Time",fontsize=15)

    plt.yticks(ylabel,np.round(ylabel*channel_spacing/1000+MIN_DISTANCE,2),size=15)
    plt.ylabel("Distance(Km)",fontsize=15)
    plt.axvline(int(Pass_time_dim), color='red')

    plt.savefig(image_filename,bbox_inches = 'tight')
    plt.close()


#%% Mean Sampling
#Mean sampling





