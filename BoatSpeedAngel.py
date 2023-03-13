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
from math import asin,degrees,sqrt,tanh,sinh,cos
from scipy import interpolate
from numpy import NaN
#输入船只两侧散波在DAS图上的斜率，单位为m/s。
ship=0#选择不同的船只
#Ais数据
#对应的MMSI
MMSI=['413260090','413208430','413471740']
V=[13.86,13.82,12.83]
Angle=[119.86,58.22,123.26]
print(V[ship],Angle[ship])


K1=[8.40,28.9,8.42] #14.61
#K1=[8.02,28.9,8.42] #14.61

#更靠近通道4000
K2=[19.2,8.4,20.15]
#K2=[28,8.4,20.15]

#%%

k1=K1[ship]
k2=K2[ship]



def f1(fai,alpha,k1,k2):
    err=k1*sin(radians(fai-alpha))-k2*sin(radians(180-fai-alpha))
    '''
    if sin(radians(180-fai-alpha))==0:
        err=1
    else:
        err=(k1*sin(radians(fai-alpha))/(k2*sin(radians(180-fai-alpha)))-1)
    '''
    return err

res=pd.DataFrame(columns=['fai','alpha','err1','speed']) #speed 指代散波的前进速度

for fai in tqdm(np.arange(1,180,0.2)):
    for a in np.arange(1,60,0.2):
        if k1==0:
            #if a<=90:
            res.loc[len(res.index)] = (fai,a,abs(a-fai),k2*sin(radians(fai+a)))
        elif k2==0:
            #if a>90:
            res.loc[len(res.index)] = (fai,a,abs(180-(a+fai)),k1*sin(radians(2*a)))                
        else:
            res.loc[len(res.index)] = (fai,a,(f1(fai,a,k1,k2)),k1*sin(radians(fai-a)))

g=9.8

#%%按照水深H与船速来筛选可能的航速和航向，考虑FR大于1

res['Height']=res['speed']**2/g
res1=res[((res['err1'])>-0.02)& ((res['err1'])<=0.02) & (res['Height']>=7)&(res['Height']<9)]




#%%考虑kelvin wedge 包络线所对应的夹角



fai=list(res1['fai'])
alpha=list(res1['alpha'])

K_CUSP=[4.6,4.5,3.67]
k_cusp=K_CUSP[ship]
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
    ShipSpeed_r.append(k_cusp*sin(radians(fai[i]-r_list[i]))/sin(radians(r_list[i])))

res1['ShipSpeed_env']=ShipSpeed_r #利用包络线计算船速
#res1['ShipSpeed(Kn)']=res1['ShipSpeed_r']/0.51444

res1['R']=r_list
res1['Fr1']=res1['ShipSpeed_env']/res1['speed'] #求解FR
res1['Fr0']=FR_list   #根据alpha-r曲线寻找Fr
res1['ShipSpeed']=res1['speed']*res1['Fr0'] #利用Fr计算船速
if ship==1:
    k=k2
else:
    k=k1

Fr1=list(res1['Fr1'])

V_error=[]
for i in range(0,len(fai)):
    V_error.append(1-k_cusp/k/(Fr1[i]*sin(radians(fai[i]))*sin(radians(fai[i]-alpha[i]))/sin(fai[i]-r_list[i])))

res1['V_error']=V_error
res1=res1[(res1['ShipSpeed_env']<20)]


res1.to_excel('船只轨迹计算结果'+'.xlsx')

#%%考虑Fr少于1的情况

