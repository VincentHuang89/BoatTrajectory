import math
from math import sin,pi,radians
from re import U
import numpy as np
from numpy import angle
import pandas as pd
from tqdm import tqdm
from WaveVAngel import FroudeNum,FRandAngel
from math import acos,degrees
import matplotlib.pyplot as plt
from math import asin,degrees,sqrt,tanh,sinh,cos
from scipy import interpolate
from numpy import NaN
import os

DAS_Results_url = '/home/huangwj/DAS/BoatTrajectory/413260090船只轨迹计算结果.xlsx'
df = pd.read_excel(DAS_Results_url, index_col=False, header=0)

fig = plt.figure(dpi=500,figsize=(8, 6))

ax1 = fig.add_subplot(111)
lns1 = ax1.plot(df['Height'], df['phi'], '.-', color='b', label=r'$\varphi$')
# ax1.plot(7.298,119.8,'y*')
ax1.set_ylabel(r'Angle ($^\circ$)')
ax1.set_xlabel('Water depth (m)',fontsize=12)

ax2 = ax1.twinx()
lns3 = ax2.plot(df['Height'], df['ShipSpeed_env'], '.-', color='r', label='ship speed')
# ax2.plot(7.298,14.02,'y*')
ax2.set_ylabel('Speed (m/s)')
ax2.set_xlabel('same')

lns = lns3 + lns1
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc=0, fontsize=14)
plt.savefig('/home/huangwj/DAS/BoatTrajectory/Paperfig/ss.pdf',bbox_inches='tight')
