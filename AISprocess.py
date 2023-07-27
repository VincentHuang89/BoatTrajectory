
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
import os
os.chdir('/home/huangwj/DAS/BoatTrajectory/')


flag=2
if flag ==1:
#2022 7月15到8月15 AIS数据
    print("处理2022 7月15到8月15 AIS数据")
    PosFile='AIS_DATA_RAW/pos_zf_gsd_1657814400_1660579199_742.csv'
    StaticFile='AIS_DATA_RAW/static_zf_gsd_1657814400_1660579199_742.csv'
    ST_UTC8=datetime.datetime.strptime("15/07/22 00:01", "%d/%m/%y %H:%M")
    ET_UTC8=datetime.datetime.strptime("15/08/22 23:59", "%d/%m/%y %H:%M")
    filename='AIS_SHIP_DATA/FiberBoatMessage_220715_220815.csv'
elif flag ==2:
#2021 9月份AIS数据
    print('处理2021 9月份AIS数据')
    PosFile='AIS_DATA_RAW/pos_zf_gsd_2_1630425600_1633017599_743.csv'
    StaticFile='AIS_DATA_RAW/static_zf_gsd_2_1630425600_1633017599_743.csv' 
    ST_UTC8=datetime.datetime.strptime("01/09/21 00:01", "%d/%m/%y %H:%M")
    ET_UTC8=datetime.datetime.strptime("30/09/21 23:59", "%d/%m/%y %H:%M")
    filename='AIS_SHIP_DATA/FiberBoatMessage_210901_210930.csv'

elif flag ==3:
#2021 10-01到12-01 的AIS数据
    print('处理2021 10-01到12-01 的AIS数据 ')
    PosFile='AIS_DATA_RAW/pos_guishan_v2_20211001120000_20211201120000_2402.csv'
    StaticFile='AIS_DATA_RAW/static_guishan_v2_20211001120000_20211201120000_2402.csv' 
    ST_UTC8=datetime.datetime.strptime("01/10/21 00:01", "%d/%m/%y %H:%M")
    ET_UTC8=datetime.datetime.strptime("01/12/21 23:59", "%d/%m/%y %H:%M")
    filename='AIS_SHIP_DATA/FiberBoatMessage_211001_211201.csv'

elif flag ==4:
#2022 04-01  07-01 的AIS数据
    print('处理2022 04-01  07-01 的AIS数据')
    PosFile='AIS_DATA_RAW/pos_guishan_v3_20220401120000_20220701120000_2403.csv'
    StaticFile='AIS_DATA_RAW/static_guishan_v3_20220401120000_20220701120000_2403.csv' 
    ST_UTC8=datetime.datetime.strptime("01/04/22 00:01", "%d/%m/%y %H:%M")
    ET_UTC8=datetime.datetime.strptime("01/07/22 23:59", "%d/%m/%y %H:%M")
    filename='AIS_SHIP_DATA/FiberBoatMessage_220401_220701.csv'

else:
    print('处理22年9月到23年5月的数据')
    PosFile='AIS_DATA_RAW/pos_guishan_20220901120000_20230501120000_2401.csv'
    StaticFile='AIS_DATA_RAW/static_guishan_20220901120000_20230501120000_2401.csv' 
    ST_UTC8=datetime.datetime.strptime("01/09/22 00:01", "%d/%m/%y %H:%M")
    ET_UTC8=datetime.datetime.strptime("01/05/23 23:59", "%d/%m/%y %H:%M")
    filename='AIS_SHIP_DATA/FiberBoatMessage_220901_230501.csv'




df=pd.read_csv(StaticFile)
dfp=pd.read_csv(PosFile)
print(f'statis:{len(df)}',f'pos:{len(dfp)}')
print(df.columns)
print(dfp.columns)

ST_UTC0=ST_UTC8-timedelta(hours=8)
ET_UTC0=ET_UTC8-timedelta(hours=8)

FiberBoatMessage=AISData(PosFile,StaticFile,ST_UTC8,ET_UTC8)
FiberBoatMessage.sort_values(by='CrossTime',ascending=True,inplace=True,ignore_index=True)
FiberBoatMessage.to_csv(filename)
print(f'Write {filename}')

#%%画出船轨迹
'''
M= folium.Map(location=[22.1,113.8],width=600,height=600,zoom_start=12)
MMSI=['413260090'] 
PlotTraj(PosFile,StaticFile,ST_UTC8,ET_UTC8,MMSI,M)

#%%
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
'''