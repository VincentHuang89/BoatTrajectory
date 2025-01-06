#%%依照时间约束在OSM中显示过光纤的船只轨迹
import pandas as pd
from IPython.display import display
import matplotlib.pyplot as plt
#import coordTransform as ct
import numpy as np
from scipy import stats
import pyproj 
from sklearn.linear_model import RANSACRegressor
import os
import math
import datetime
from datetime import timedelta
from math import sin,cos
import folium
import os
from tqdm import tqdm
from shapely.geometry import LineString


def convert_to_cartesian(latitude, longitude, radius= 6378.137):
    """
    将经纬度坐标转换为直角坐标系（笛卡尔坐标系）
    
    参数：
    latitude: 纬度，单位为度数
    longitude: 经度，单位为度数
    radius: 地球半径
    
    返回值：
    包含 x、y、z 直角坐标的 numpy 数组
    """
    lat_rad = np.radians(latitude)
    lon_rad = np.radians(longitude)
    
    x = radius * np.cos(lat_rad) * np.cos(lon_rad)
    y = radius * np.cos(lat_rad) * np.sin(lon_rad)
    z = radius * np.sin(lat_rad)
    
    return np.array([x, y, z])


def cartesian_to_geo(x, y, z):
    """
    Convert cartesian coordinates to geographical coordinates
    
    Args:
        x: ECEF x coordinate 
        y: ECEF y coordinate
        z: ECEF z coordinate
        
    Returns:
        (lon, lat): Tuple of longitude and latitude
    """
    # Calculate horizontal distance
    xy_dist = math.sqrt(x**2 + y**2)
    # Calculate latitude
    lat = math.atan2(z, xy_dist) 
    # Calculate longitude 
    lon = math.atan2(y, x)
    # Convert from radians to degrees
    lat = math.degrees(lat)
    lon = math.degrees(lon)
    
    return lon, lat
def cross(A,B,C,D):
    track1=[A,B]
    track2=[C,D]
    line1 = LineString(track1)
    line2 = LineString(track2)
    
    if line1.intersects(line2):
        return True
    else:
        return False

FIBER_0=(22.161561, 113.712681)  #纬度（latitude）  经度  （longtitude） （调整后光纤位置）
FIBER_1=(22.167286, 113.795594) #桂山岛方向 

filepath = '/home/huangwj/DAS/BoatTrajectory/'

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


df=pd.read_csv(os.path.join(filepath,StaticFile))
dfp=pd.read_csv(os.path.join(filepath,PosFile))
dfp.rename(columns={'MMSI':'MMSI','经度':'lng','纬度':'lat','航向':'heading','航速(节)':'speed','航艏向':'Ship heading','时间':'Time'},inplace=True)

dfp = dfp[dfp['Time'].notnull()]
print(dfp.dtypes)
dfp['Time'] = pd.to_datetime(dfp['Time'],errors='ignore',format = '%Y-%m-%d %H:%M:%S')
dfp_0 = dfp[dfp['MMSI'] == 413696170]
print(f'statis:{len(df)}',f'pos:{len(dfp)}')
#print(df.columns)
print(dfp.columns)

ST_UTC0=ST_UTC8-timedelta(hours=8)
ET_UTC0=ET_UTC8-timedelta(hours=8)
print(dfp.dtypes)


#依照时间筛选船只的AIS数据
Time_start = '01/09/21 00:01'
Time_end = '01/09/21 00:30'
ts =datetime.datetime.strptime(Time_start, "%d/%m/%y %H:%M")
te =datetime.datetime.strptime(Time_end, "%d/%m/%y %H:%M")
dfp = dfp[dfp['Time']>=ts]
dfp = dfp[dfp['Time']<=te]
#%%
print(set(dfp['MMSI'].values))
Ship_set = set(dfp['MMSI'].values)
#筛选过光纤的船轨迹
Ship_cross_fiber=[]
for ship in Ship_set:
    dfp_ship = dfp[dfp['MMSI']==ship]
    dfp_ship.sort_values(by=['Time'],ascending = True, inplace = True)
    lng_list = dfp_ship['lng'].tolist()
    lat_list = dfp_ship['lat'].tolist()
    for i in range(0,len(lng_list)-1):
        Tra_0=np.array([lat_list[i],lng_list[i]])
        Tra_1=np.array([lat_list[i+1],lng_list[i+1]])
        if cross(FIBER_0,FIBER_1,Tra_0,Tra_1)==1:
            Ship_cross_fiber.append(ship)


Ship_cross_fiber = set(Ship_cross_fiber)

m = folium.Map(location=[22.1,113.8],width=600,height=600,zoom_start=12)
folium.PolyLine([[FIBER_0[0],FIBER_0[1]],[FIBER_1[0],FIBER_1[1]]],color='green').add_to(m)
for ship in Ship_cross_fiber:
    dfp_ship = dfp[dfp['MMSI']==ship]
    dfp_ship.sort_values(by=['Time'],ascending = True, inplace = True)
    lng_list = dfp_ship['lng'].tolist()
    lat_list = dfp_ship['lat'].tolist()
    for i in range(0,len(lng_list)-1):
        folium.PolyLine([[lat_list[i],lng_list[i]],[lat_list[i+1],lng_list[i+1]]]).add_to(m)

display(m)
m.save('map.html')


#%%

#读取csv
os.chdir('/home/huangwj/DAS/DataSet/')
POS_TIME = pd.read_csv('merge_pos_time.csv',index_col=0)
POS_TIME['t0'] = pd.to_datetime(POS_TIME['t0'],errors='ignore',format = '%Y-%m-%d %H:%M:%S')
POS_TIME['t1'] = pd.to_datetime(POS_TIME['t1'],errors='ignore',format = '%Y-%m-%d %H:%M:%S')
POS_TIME['Delta_AIS_TIME'] = POS_TIME['t1']-POS_TIME['t0']

#筛选时间间距较少的记录
POS_TIME['Delta_AIS_TIME'] = POS_TIME['Delta_AIS_TIME'].apply(lambda x: x.total_seconds())
POS_TIME = POS_TIME[POS_TIME['Delta_AIS_TIME']<=10]
POS_cal = POS_TIME['POS'][POS_TIME['Speed'].astype(float)>8].tolist()

#筛选不正常的数据点

POS_cal = [pos for pos in POS_cal if eval(pos)[1]<22.169]





