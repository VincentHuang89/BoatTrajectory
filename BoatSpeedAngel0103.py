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
MMSI=['413471740','413260090','413208430']
V=[12.83,13.86,13.82,7.078,13.11]
Angle=[123.26,119.86,58.22,90.707,120.08]
print(V[ship],Angle[ship])

L1=np.array([6,6,8,14,4])  #实验部分船只的数据Km
t1=np.array([11,11.4,4.6,24.5,8.25])
K1=L1*1000/(t1*60)
K1=[8.796,8.48,28.9] 
#更靠近通道4000
L2=np.array([8,8,8,16,0])
t2=np.array([4.2,5,12.8,19.2,1])
K2=L2*1000/(t2*60)
K2=[18.7,19.2,10.11]

#%%

k1=K1[ship]
k2=K2[ship]



def f1(fai,alpha,k1,k2):
    #err=k1*sin(radians(fai-alpha))-k2*sin(radians(180-fai-alpha))
    if sin(radians(180-fai-alpha))==0:
        err=1
    else:
        err=(k1*sin(radians(fai-alpha))/(k2*sin(radians(180-fai-alpha)))-1)
    return err

res=pd.DataFrame(columns=['fai','alpha','err1','speed']) #speed 指代散波的前进速度

for fai in tqdm(np.arange(1,180,0.5)):
    for a in np.arange(1,60,0.5):
        if k1==0:
            #if a<=90:
            res.loc[len(res.index)] = (fai,a,abs(a-fai),k2*sin(radians(fai+a)))
        elif k2==0:
            #if a>90:
            res.loc[len(res.index)] = (fai,a,abs(180-(a+fai)),k1*sin(radians(2*a)))                
        else:
            res.loc[len(res.index)] = (fai,a,abs(f1(fai,a,k1,k2)),k1*sin(radians(fai-a)))

g=9.8

#%%按照水深H与船速来筛选可能的航速和航向，考虑FR大于1

res['Height']=res['speed']**2/g
res1=res[(abs(res['err1'])<0.01) & (res['Height']>=7)]




#%%考虑kelvin wedge 包络线所对应的夹角



fai=list(res1['fai'])
alpha=list(res1['alpha'])
speed=list(res1['speed'])
k_cusp=3.67
k_cusp=4.6
#k_cusp=4.5
#fai=list(180-np.array(fai))

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

#根据res1中的alpha，结合Alpha-r曲线来修订平行等距波峰线的alpha

def SearchParam(alpha,DF):
    '''
    返回DF中大于输入alpha的所有alpha，R值和Fr值
    '''
    Alpha=list(DF['Alpha'])
    R=list(DF['R'])
    Fr=list(DF['FroudeNumber'])
    alpha_err=list(np.abs(np.array(Alpha)-alpha))
    pos=alpha_err.index(min(alpha_err))
    return Alpha[0:pos],R[0:pos],Fr[0:pos]

def FindAlpha(alpha_range,R_range,Fr_range,speed,k_cusp,fai):
    '''
     在alpha_range中寻找合适的alpha
    '''
    r_list=np.array(R_range)
    s1=k_cusp*np.sin(np.radians(fai-r_list))/np.sin(np.radians(r_list))

    s2=speed*np.array(Fr_range)    
    print(s1)
    print(s2)
    err=list(s1-s2)
    pos=err.index(min(err))
    return alpha_range[pos],R_range[pos],Fr_range[pos]



r_list=[]
FR_list=[]
alpha_list=[]
for i in range(0,len(alpha)):
    alpha_range,R_range,Fr_range=SearchParam(alpha[i],DF)
    alpha_rev,r_rev,Fr_rev=FindAlpha(alpha_range,R_range,Fr_range,speed[i],k_cusp,fai[i])
    r_list.append(r_rev)
    FR_list.append(Fr_rev)
    alpha_list.append(alpha_rev)


ShipSpeed_r=[]
for i in range(0,len(fai)):
    ShipSpeed_r.append(k_cusp*sin(radians(fai[i]-r_list[i]))/sin(radians(r_list[i])))

res1['ShipSpeed_r']=ShipSpeed_r #包络线计算速度
#res1['ShipSpeed(Kn)']=res1['ShipSpeed_r']/0.51444

res1['R']=r_list
res1['Fr1']=res1['ShipSpeed_r']/res1['speed']
res1['Fr0']=FR_list
res1['alpha0']=alpha_list
res1['ShipSpeed']=res1['speed']*res1['Fr0'] #利用apha索引Fr，计算速度
res1['ShipSpeed_0']=np.array(res1['speed'])/np.sin(np.radians(res1['alpha'])) #利用V/sin(alpha)
res1.to_excel('船只轨迹计算结果'+'.xlsx')

#%%考虑Fr少于1的情况

