#在时间维度上对齐AIS数据与DAS实测数据

#增加datetime类型
#%%
import numpy as np
import matplotlib.pyplot as plt
from numpy import NaN
import pandas as pd
from obspy import UTCDateTime

df_JT = pd.read_excel('JIANTING.xlsx')
df_JT['time']= df_JT['time'].astype(str)

print(df_JT)
#%%
def  ExtractTIme(TimeStr:str):
    Year=[]
    Mon=[]
    Day=[]
    Hour=[]
    Min=[]
    Sec=[]
    for i in range(0,len(TimeStr)):
        Year.append(int(TimeStr[i][0:4]))
        Mon.append(int(TimeStr[i][5:7]))
        Day.append(int(TimeStr[i][8:10]))
        Hour.append(int(TimeStr[i][11:13]))
        Min.append(int(TimeStr[i][14:16]))
        Sec.append(int(TimeStr[i][17:19]))
    return Year,Mon,Day,Hour,Min,Sec


Year,Mon,Day,Hour,Min,Sec = ExtractTIme(df_JT['time'])
print(Year,Mon,Day,Hour,Min,Sec)
# %%计算每一秒钟的经纬度信息
#temp=np.transpose((np.concatenate([Year,Mon,Day,Hour,Min,Sec],axis=0)).reshape(6,14))
temp=np.transpose((np.concatenate([Hour,Min,Sec],axis=0)).reshape(3,14))
#%%
#DF_JT=pd.DataFrame(temp,columns=['Year','Mon','Day','Hour','Min','Sec'])
DF_JT=pd.DataFrame(temp,columns=['Hour','Min','Sec'])
#DF_JT.insert(DF_JT.shape[1],'Year',Year)
#DF_JT=pd.DataFrame([Year,Mon,Day,Hour,Min,Sec],columns=['Year','Mon','Day','Hour','Min','Sec'])
DF_JT=pd.concat([df_JT[['longtitude','latitude']],DF_JT],axis=1)

print(DF_JT)

# %%
DF_JT_Quan= pd.DataFrame([],columns=DF_JT.columns)
print(DF_JT_Quan)
print(DF_JT.shape[0])

#%%
for i in range(0,DF_JT.shape[0]-1):
    DF_JT_Quan=DF_JT_Quan.append(DF_JT.iloc[i])
    Hour_diff=DF_JT.iloc[i+1]['Hour']-DF_JT.iloc[i]['Hour']
    Min_diff=DF_JT.iloc[i+1]['Min']-DF_JT.iloc[i]['Min']
    Sec_diff=DF_JT.iloc[i+1]['Sec']-DF_JT.iloc[i]['Sec']
    long_diff=DF_JT.iloc[i+1]['longtitude']-DF_JT.iloc[i]['longtitude']
    lat_diff=DF_JT.iloc[i+1]['latitude']-DF_JT.iloc[i]['latitude']
    Time_diff= int(Hour_diff*60+Min_diff*60+Sec_diff)
    for j in range(0,Time_diff-1):
        long=DF_JT.iloc[i]['longtitude']+(j+1)*long_diff/Time_diff
        lat=DF_JT.iloc[i]['latitude']+(j+1)*lat_diff/Time_diff
        DF_JT_Quan= DF_JT_Quan.append({'longtitude':long,'latitude':lat},ignore_index=True)
DF_JT_Quan=DF_JT_Quan.append(DF_JT.iloc[DF_JT.shape[0]-1],ignore_index=True)

print(DF_JT_Quan)


# %%
