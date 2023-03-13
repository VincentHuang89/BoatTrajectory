from turtle import color
import numpy as np
import pandas as pd
import math
import datetime
from datetime import timedelta
from math import sin,cos
import folium


def PlotTraj(PosFile:str,StaticFile:str,ST_UTC8:datetime.datetime,ET_UTC8:datetime.datetime,MMSI:list,M):
    df = pd.read_csv(PosFile) #UTC+8
    static=pd.read_csv(StaticFile)

    df['MMSI'].value_counts()
    df['时间']=pd.to_datetime(df['时间'])
    df_time=pd.DataFrame(df[(df['时间']>ST_UTC8) &  (df['时间']<ET_UTC8)])
    INDEX=df_time['MMSI'].value_counts().index
    temp=df_time['MMSI'].value_counts()

    #删除只有一个数据的轨迹
    dropindex=[]
    for i in range(0,len(temp)):
        if temp.loc[INDEX[i]]==1:
            dropindex.append(temp.index[i])
    df_time=df_time[~df_time['MMSI'].isin(dropindex)]
    df_time['MMSI'] = df_time['MMSI'].astype(str)
    

    colobars=['black','yellow','blue','green','white']
    i=0
    for mmsi in MMSI:
        df_tra=df_time[df_time['MMSI']==mmsi]
        tra=np.transpose([df_tra['纬度'].values,df_tra['经度'].values])
        folium.PolyLine(tra,color=colobars[i]).add_to(M)
        i=i+1
    folium.PolyLine([[22.140,113.709],[22.168,113.801]],color='red').add_to(M)
    return M
