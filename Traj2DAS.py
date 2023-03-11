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
phy=phy-90
dist=15

N=FroudeNum(v,h)

X,Y,Y1,alpha=WavePattern1(1,N)
PlotShipWakePattern(v,h,15,phy,dist)


PlotBoatWave(v,h,phy,dist)

#%%
t=5

PlotSimuWaveInDas(v,h,phy,t,8,1,0)
#%%
WLen_Scale=35
Wbias=73.5
Tbias=1.2
h=7.125
angle=118.2

v=13.262
UpperBound=20
LowerBound=-1000

PlotSimuWaveInDas(v,h,angle,t,12,WLen_Scale,Wbias)
# %%

