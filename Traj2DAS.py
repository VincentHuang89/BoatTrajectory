#%%
'''
Simulate the divergent wave perceived by DAS according to the classic ship wake theory
'''
from cmath import isnan
from matplotlib.lines import lineStyles
import matplotlib.pyplot as plt
import numpy as np
import math
from math import sin,pi,radians,tan,cos,sin
import numpy as np
import pandas as pd
from WaveVAngel import FroudeNum,FRandAngel,WavePattern,ROTATE,move,cross,WaveCross,PlotWaveInDAS,PlotWaveInDas,PlotBoatWave,PlotSimuWaveInDas,WavePattern1,PlotShipWakePattern
from math import acos,degrees
from PIL import Image

#%% 
v=13.02
h=7.298
phy=120   #船
dist=10

N=FroudeNum(v,h)
#X,Y,Y1,alpha=WavePattern1(1,N)  
#PlotShipWakePattern(v,h,5,phy,dist)

PlotBoatWave(v,h,phy,dist)  #画出不同时刻船只经过光纤的过程

#%%
WLen_Scale=1
Wbias=0
Tbias=0
angle=phy
t=10

UpperBound=40
LowerBound=-1000

PlotSimuWaveInDas(v,h,angle,t,8,WLen_Scale,Wbias,UpperBound,LowerBound)
# %%

