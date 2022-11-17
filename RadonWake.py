
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
from math import sin,cos,tan,radians
from sklearn import preprocessing
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

import skimage 
from skimage.transform import radon
def PlotRadon(ShowDataSlice):
    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1.set_title("Original")
    ax1.imshow(ShowDataSlice, aspect='auto',cmap="bwr",origin='lower',vmin=-3,vmax=3)
    theta = np.linspace(0., 180., max(ShowDataSlice.shape), endpoint=False)
    sinogram = radon(ShowDataSlice, theta=theta)
    sinogram=(sinogram-np.mean(sinogram))/np.std(sinogram,ddof=1)
    STD=3*np.std(sinogram,ddof=1)
    sinogram=((sinogram>STD)|(sinogram<-STD))*sinogram
    dx, dy = 0.5 * 180.0 / max(ShowDataSlice.shape), 0.5 / sinogram.shape[0]
    ax2.set_title("Radon transform\n(Sinogram)")
    ax2.set_xlabel("Projection angle (deg)")
    ax2.set_ylabel("Projection position (pixels)")
    SinSlice=slice(int(0.45*sinogram.shape[1]),int(0.55*sinogram.shape[1]),1)
    #ax2.imshow(sinogram[:,SinSlice], cmap="bwr",aspect='auto',origin='lower')
    ax2.imshow(sinogram, cmap="bwr",extent=(-dx, 180.0 + dx, -dy, sinogram.shape[0] + dy),aspect='auto',origin='lower')
    plt.savefig('Radon_transform.png')
    return sinogram

def SpeedOnRadon(sinogram,resolution,channel_spacing,fs):
#求解最大的斜率
    deg_y=np.argmax(sinogram)%sinogram.shape[1]/resolution*180
    pos_x=int(np.argmax(sinogram)/sinogram.shape[1])
    print(pos_x,deg_y)
    #将斜率换算成速度m/s
    speed=1/tan(radians(deg_y))*channel_spacing*fs
    #speed1=tan(radians(90-deg_y))*channel_spacing*fs
    print('Radon 域最大值所对应的deg相应的速度',speed)

    #求解最大的斜率，平均每行的最大值索引
    deg_y_list=[]
    sino=[]
    for i in range(0,sinogram.shape[0]):
        sino.append(np.max(sinogram[i,:]))
        deg_y_list.append(np.argmax(sinogram[i,:])%sinogram.shape[1]/resolution*180)
    df=pd.DataFrame(np.transpose(np.array([sino,deg_y_list])),columns=['sino','deg'])
    df['sino']=abs(df['sino'])
    df.sort_values(by='sino',ascending=False,inplace=True)
    deg_y_list=df['deg'][0:10]
    print(np.mean(deg_y_list),np.median(deg_y_list))
    s1=1/tan(radians(np.mean(deg_y_list)))*channel_spacing*fs
    s2=1/tan(radians(np.median(deg_y_list)))*channel_spacing*fs
    print('选取Radon域上前10的最大值所对应的deg\n')
    print('均值',s1,'中位数',s2)

def SNROnRadon(sinogram):
    #[1]M. T. Rey, J. K. E. Tunaley, J. T. Folinsbee, P. A. Jahans, J. A. Dixon, and M. R. Vant, "Application Of Radon Transform Techniques To Wake Detection In Seasat-A SAR Images," IEEE Transactions on Geoscience and Remote Sensing, vol. 28, pp. 553-560, 1990.

    print('SNR on Radon domain: ',(np.max(sinogram)-np.mean(sinogram))/np.std(sinogram,ddof=1))
