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

#输入船只两侧散波在DAS图上的斜率，单位为m/s，是相速度，不是群速度
ship=1
#Ais数据
V=[12.83,13.86,13.82,7.078,13.11]
Angle=[123.26,119.86,58.22,90.707,120.08]
print(V[ship],Angle[ship])


L1=np.array([6,6,8,14,4])  #Km
t1=np.array([11,11.4,4.6,24.5,8.25])
K1=L1*1000/(t1*60)

#更靠近通道4000
L2=np.array([8,8,8,16,0])
t2=np.array([4.2,5,12.8,19.2,1])
K2=L2*1000/(t2*60)

#%%

k1=K1[ship]
k2=K2[ship]

def f1(fai,alpha,k1,k2):
    err=k1*sin(radians(fai-alpha))-k2*sin(radians(180-fai-alpha))
    return err



res=pd.DataFrame(columns=['fai','alpha','err1','speed'])

for fai in tqdm(np.arange(1,180,1)):
    for a in np.arange(0,90,1):
        if k1==0:
            #if a<=90:
            res.loc[len(res.index)] = (fai,a,abs(a-fai),k2*sin(radians(fai+a)))
        elif k2==0:
            #if a>90:
            res.loc[len(res.index)] = (fai,a,abs(180-(a+fai)),k1*sin(radians(2*a)))                
        else:
            res.loc[len(res.index)] = (fai,a,abs(f1(fai,a,k1,k2)),k1*sin(radians(fai-a)))

g=9.8

#%%

res['Height']=res['speed']**2/g
res1=res[(abs(res['err1'])<0.02) & (res['Height']>=7)]
res1['ShipSpeed']=res1['speed']/np.sin(np.deg2rad(res1['alpha']))
res1['ShipSpeed(Kn)']=res1['ShipSpeed']/0.51444
res1['Fr']=res1['ShipSpeed']/res1['speed']

#res1['BoatSpeed']=res1['speed']/np.sin(np.radians(np.array(res1['beta'])))

#%%
res1.to_excel('船只轨迹计算结果'+str(U[ship])+'-'+str(Angle[ship])+'.xlsx')

#%%分析船只的佛汝德系数范围
print(
FroudeNum(6,7),
FroudeNum(6,14),
FroudeNum(12,7),
FroudeNum(12,10))

#%%佛汝德系数对alpha角度的影响
Fr=FroudeNum(12,10)
print(degrees(asin(1/Fr)))
Fr=FroudeNum(12,7)
print(degrees(asin(1/Fr)))
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
for Fr in np.arange(1,3,0.01):
        FR.append(Fr)
        alpha=degrees(asin(1/Fr))
        ALPHA.append(alpha)
plt.figure(dpi=300)
plt.plot(FR,ALPHA)
plt.xlabel('Fr')
plt.ylabel(chr(945))

# %%
degrees(asin(1/3))
# %%

# %%假设单边波峰线平行于光纤
L1=6  #Km
t1=11 #Min  #25.5,26.5,11.5,25.6
k1=L1*1000/(t1*60)
Thai=123
V=12.82#13.856#12.83
H=7.4
g=9.8
fr=FroudeNum(V,H)
alpha1=180-Thai-degrees(asin(sqrt(g*H)/k1))
alpha1=Thai-degrees(asin(sqrt(g*H)/k1))
alpha2=degrees(asin(1/sqrt(fr)))

print(alpha1,alpha2)
#%%
V=10
H=np.arange(6,8,0.01)
FR=[]
Alpha=[]
for h in H:
    fr=FroudeNum(V,h)
    FR.append(fr)
    alpha2=degrees(asin(1/sqrt(fr)))
    Alpha.append(alpha2)
plt.figure()
plt.plot(FR,Alpha)
plt.xlabel('Fr')
plt.ylabel('Degree')
# %%
V=13
h=7
fr=FroudeNum(V,h)



# %%
