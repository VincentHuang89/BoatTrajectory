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
from math import asin,degrees,sqrt,tanh,sinh
from scipy import interpolate
from numpy import NaN
#输入船只两侧散波在DAS图上的斜率，单位为m/s。
ship=0  #选择不同的船只
#Ais数据
#对应的MMSI
MMSI=['413471740','413260090','413208430']
V=[12.83,13.86,13.82,7.078,13.11]
Angle=[123.26,119.86,58.22,90.707,120.08]
print(V[ship],Angle[ship])

L1=np.array([6,6,8,14,4])  #实验部分船只的数据Km
t1=np.array([11,11.4,4.6,24.5,8.25])
K1=L1*1000/(t1*60)
K1=[8.476,8.304,14.97]

#更靠近通道4000
L2=np.array([8,8,8,16,0])
t2=np.array([4.2,5,12.8,19.2,1])
K2=L2*1000/(t2*60)
K2=[18,18.512,9.406]
#%%

k1=K1[ship]
k2=K2[ship]



def f1(fai,alpha,k1,k2):
    err=k1*sin(radians(fai-alpha))-k2*sin(radians(180-fai-alpha))
    return err

res=pd.DataFrame(columns=['fai','alpha','err1','speed']) #speed 指代散波的前进速度

for fai in tqdm(np.arange(1,180,1)):
    for a in np.arange(1,90,1):
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
res1=res[(abs(res['err1'])<0.02) & (res['Height']>=6)]
res1['ShipSpeed']=res1['speed']/np.sin(np.deg2rad(res1['alpha']))
res1['ShipSpeed(Kn)']=res1['ShipSpeed']/0.51444
res1['Fr']=res1['ShipSpeed']/res1['speed']
res1=res1[res1['Fr']<=10]
res1.to_excel('船只轨迹计算结果'+'.xlsx')

#%%考虑Fr少于1的情况

res2=res[(abs(res['err1'])<0.2)]
res2=res2[res2['alpha']>=10]
res2.drop(columns=['Height'],inplace=True)
res2['ShipSpeed']=res2['speed']/np.sin(np.deg2rad(res2['alpha']))
res2['ShipSpeed(Kn)']=res2['ShipSpeed']/0.51444


#res1['BoatSpeed']=res1['speed']/np.sin(np.radians(np.array(res1['beta'])))

#%%



#%%画出 curve of Fr and alpha

H=7
ALPHA=[]
FR=[]
for k in np.arange(0.01,2.5,0.01):
    n=2*k*H/sinh(2*k*H)
    V=(sqrt((tanh(k*H)/(k*H)*(3-n)*g*H)/2))
    Fr=V/(sqrt(g*H))
    FR.append(Fr)
    if V<(sqrt(g*H)):
        alpha= degrees(acos((sqrt(8*(1-n))/(3-n))))
        ALPHA.append(alpha)
FR.sort()
ALPHA.sort()
for Fr in np.arange(1,10,0.01):
        FR.append(Fr)
        alpha=degrees(asin(1/Fr))
        ALPHA.append(alpha)
plt.figure(dpi=300)
plt.plot(FR,ALPHA)
plt.xlabel('Fr')
plt.ylabel(chr(945))
#

#判断alpha与Fr的关系是否服从Fr-alpha曲线
f = interpolate.interp1d(FR, ALPHA, kind='linear')
fr=list(np.arange(np.min(FR),np.max(FR),0.01))
alpha=list(np.round(f(fr),4))
plt.plot(fr,alpha)
#%%
fr=list(np.round(fr,2))
alpha[fr.index(1.21)]
FR_find=list(np.round(res1['Fr'],2))
a=[]
for i in range(0,len(FR_find)):
    a.append(alpha[fr.index(FR_find[i])])
res1['alpha1']=a
# %%重新插值，将alpha作为X轴，Fr作为Y轴
AlphaLess1=[]
FrLess1=[]
for i in range(0,len(FR)):
    if FR[i]<=1:
        AlphaLess1.append(ALPHA[i])
        FrLess1.append(FR[i])
AlphaMore1=[]
FrMore1=[]
for i in range(0,len(FR)):
    if FR[i]>1:
        AlphaMore1.append(ALPHA[i])
        FrMore1.append(FR[i])

f0 = interpolate.interp1d(AlphaLess1, FrLess1, kind='linear')
alpha0=list(np.arange(np.min(AlphaLess1),np.max(AlphaLess1),0.0001))
fr0=list(np.round(f0(alpha0),4))
alpha0=list(np.round(alpha0,4))


f1 = interpolate.interp1d(AlphaMore1, FrMore1, kind='linear')
alpha1=list(np.arange(np.min(AlphaMore1),np.max(AlphaMore1),0.0001))
fr1=list(np.round(f1(alpha1),4))
alpha1=list(np.round(alpha1,4))



#%%
a_find=list(np.round(res2['alpha'],4))
fr_find=[]
for i in range(0,len(a_find)):
    if a_find[i]>=min(alpha0) and a_find[i]<=max(alpha0):
        fr_find.append(fr0[alpha0.index(a_find[i])])
    else:
        fr_find.append(NaN)
res2['Fr0']=fr_find

fr_find=[]
for i in range(0,len(a_find)):
    if a_find[i]>=min(alpha1) and a_find[i]<=max(alpha1):
        fr_find.append(fr1[alpha1.index(a_find[i])])
    else:
        fr_find.append(NaN)
res2['Fr1']=fr_find

#%%

# %%
fr0=np.array(res2['Fr0'])
height0=[]
shipspeed=np.array(res2['ShipSpeed'])
for i in range(0,len(fr0)):
    if np.isnan(fr0[i]):
        height0.append(NaN)
    else:
        height0.append((shipspeed[i]/fr0[i])**2/g)
res2['height0']=height0


fr1=np.array(res2['Fr1'])
height1=[]
shipspeed=np.array(res2['ShipSpeed'])
for i in range(0,len(fr1)):
    if np.isnan(fr1[i]):
        height1.append(NaN)
    else:
        height1.append((shipspeed[i]/fr1[i])**2/g)
res2['height1']=height1

# %%
