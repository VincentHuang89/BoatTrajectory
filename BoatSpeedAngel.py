#%%计算船只过光纤的夹角与速度
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
from math import asin,degrees,sqrt,tanh,sinh,cos,sin, radians
from scipy import interpolate
from numpy import NaN
#输入船只两侧散波在DAS图上的斜率，单位为m/s。
ship=3  #选择不同的船只
#Ais数据
#对应的MMSI
MMSI=['413260090','413208430','413471740','413231470','413260090']
'''V=[13.86,13.82,12.83,16.15]
Angle=[119.86,58.22,123.26,90.01]'''

V=[13.89,13.85,12.86,16.15,14.33]
Angle=[136.58,47.1,133.07,90.05,49]

print(V[ship],Angle[ship])

K1=[8.40,28.9,8.42,10.21,22.5] #8.40,
#更靠近通道4000
K2=[19.7,8.4,20.15,12.1,9.3]  #11.25


K1=[8.40,28.9,8.42,9.76,22.5] #8.40,
#更靠近通道4000
K2=[19.7,8.4,20.15,11.57,8.6]  #11.25


K_CUSP=[4.6,4.5,3.67,3.3,4.36]#3.48
k_cusp=K_CUSP[ship]


k1 = K1[ship]  # 请根据实际情况设置K1[ship]的值
k2 = K2[ship]  # 请根据实际情况设置K2[ship]的值

res = pd.DataFrame(columns=['fai', 'alpha', 'speed'])  # speed 指代散波的前进速度

alpha = np.arange(20, 60, 0.1)

K = (k1-k2)/(k1+k2)
fai = np.degrees(np.arctan(np.tan(np.radians(alpha))/K))
fai = (fai + 180) % 180
speed = k2 * np.sin(np.radians(fai + alpha))
#%%
res = pd.DataFrame({'fai': fai.flatten(), 'alpha': alpha.flatten(), 'speed': speed.flatten()})

#%%按照水深H与船速来筛选可能的航速和航向，考虑FR大于1
g=9.8
res['Height']=res['speed']**2/g
res1=res[(res['Height']>=7)&(res['Height']<10)]
#%%考虑kelvin wedge 包络线所对应的夹角

#%%

fai=list(res1['fai'])
alpha=list(res1['alpha'])


if ship==1:
    fai=list(180-np.array(fai))

#读取alpha-r曲线
DF=pd.read_excel('Alpha_R.xlsx')

#根据res1的alpha寻找合适的r
def FindR(alpha,DF):
    Alpha=list(DF['Alpha'])
    R=list(DF['R'])
    alpha_err=list(np.abs(np.array(Alpha)-alpha))
    pos=alpha_err.index(min(alpha_err))
    return(R[pos])

def FindFr(alpha,DF):
    Alpha=list(DF['Alpha'])
    Fr=list(DF['FroudeNumber'])
    alpha_err=list(np.abs(np.array(Alpha)-alpha))
    pos=alpha_err.index(min(alpha_err))
    return(Fr[pos])    

r_list=[]
FR_list=[]
for i in range(0,len(alpha)):
    r_list.append(FindR(alpha[i],DF))
    FR_list.append(FindFr(alpha[i],DF))



ShipSpeed_r=[]
for i in range(0,len(fai)):
    ShipSpeed_r.append(k_cusp*sin(radians(fai[i]-r_list[i]))/sin(radians(r_list[i])))   #考虑K=1.5时候，包络线在上半部分还是下半部分？

res1.loc[:,'ShipSpeed_env']=ShipSpeed_r #利用包络线计算船速
#res1['ShipSpeed(Kn)']=res1['ShipSpeed_r']/0.51444

res1.loc[:,'R']=r_list
res1.loc[:,'Fr1']=res1['ShipSpeed_env']/res1['speed'] #求解FR
res1.loc[:,'Fr0']=FR_list   #根据alpha-r曲线寻找Fr
res1.loc[:,'ShipSpeed']=res1['speed']*res1['Fr0'] #利用Fr计算船速
if ship==1:
    k=k2
else:
    k=k1

Fr1=list(res1['Fr1'])

V_error=[]
for i in range(0,len(fai)):
    V_error.append(1-k_cusp/k/(Fr1[i]*sin(radians(fai[i]))*sin(radians(fai[i]-alpha[i]))/sin(fai[i]-r_list[i])))

res1.loc[:,'V_error']=V_error
res1=res1[(res1['ShipSpeed_env']<20)]

print(min(list(res1['fai'])),max(list(res1['fai'])))
print(min(list(res1['ShipSpeed_env'])),max(list(res1['ShipSpeed_env'])))
print(min(list(res1['ShipSpeed'])),max(list(res1['ShipSpeed'])))

res1.to_excel('船只轨迹计算结果'+'.xlsx')
#%%
#画出解空间
Speed=list(res1['ShipSpeed_env'])
Angel=list(res1['fai'])
Depth=list(res1['Height'])

fig=plt.figure(dpi=500,figsize=(10,8))
ax1=fig.add_subplot(111)
ax1.plot(Depth,Speed,'r-',label='ship speed')
ax1.set_xlabel('Water depth (m)',fontsize=12)

ax1.set_ylabel('Speed (m/s)',fontsize=10)

ax2=ax1.twinx()
ax2.plot(Depth,Angel,'b-', label=r'\Phi')
ax2.set_ylabel('Angle ( )',fontsize=10)
plt.savefig('Paperfig/ss.pdf',bbox_inches='tight')


# %%
