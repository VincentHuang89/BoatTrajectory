
# %%
import numpy as np
import pandas as pd
import math
import datetime
from datetime import timedelta
from AISData import AISData,FilterBoat,ShipTraj
from Maptraj import PlotTraj
from math import sin,cos
import folium

flag=1

if flag ==1:
#2022 7月15到8月15 AIS数据
    PosFile='pos_zf_gsd_1657814400_1660579199_742.csv'
    StaticFile='static_zf_gsd_1657814400_1660579199_742.csv'
else:
#2021 9月份AIS数据
    PosFile='pos_zf_gsd_2_1630425600_1633017599_743.csv'
    StaticFile='static_zf_gsd_2_1630425600_1633017599_743.csv'    


#%%按照时间区间筛选
ST_UTC8=datetime.datetime.strptime("01/07/22 00:01", "%d/%m/%y %H:%M")
ET_UTC8=datetime.datetime.strptime("30/08/22 23:59", "%d/%m/%y %H:%M")
ST_UTC0=ST_UTC8-timedelta(hours=8)
ET_UTC0=ET_UTC8-timedelta(hours=8)

FiberBoatMessage=AISData(PosFile,StaticFile,ST_UTC8,ET_UTC8)
print('Data Read!')
FiberBoatMessage.sort_values(by='CrossTime',ascending=True,inplace=True,ignore_index=True)
FiberBoatMessage.to_excel('FiberBoatMessage.xlsx')

#%%画出船轨迹
'''
M= folium.Map(location=[22.1,113.8],width=600,height=600,zoom_start=12)
MMSI=['413510040','538006097'] 
M=PlotTraj(PosFile,StaticFile,ST_UTC8,ET_UTC8,MMSI,M)
M
'''
DT=timedelta(minutes=10)
Angle=[85,95]
Dist=[0,9]
Speed=[7,13]
Boat,L=FilterBoat(FiberBoatMessage,DT,Angle,Dist,Speed)
print(Boat)
Boat.to_excel('FilterBoat.xlsx')
# %%特定时间与特定船只的轨迹
ST=datetime.datetime.strptime("24/07/22 00:01", "%d/%m/%y %H:%M")
ET=datetime.datetime.strptime("24/07/22 23:59", "%d/%m/%y %H:%M")
mmsi=413260090
Df_Ship=ShipTraj(PosFile,StaticFile,ST,ET,mmsi)
Df_Ship.rename(columns={'经度':'longtitude','纬度':'latitude','航速(节)':'speed(Kn)','时间':'time'},inplace=True)
Df_Ship.drop(['航向','航艏向'],axis=1).to_excel(str(mmsi)+'.xlsx')
# %%
