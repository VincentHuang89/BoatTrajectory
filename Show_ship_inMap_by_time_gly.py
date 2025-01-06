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

    #在海图中手动拾取的经纬度信息
FIBER_0 = (24.454581,118.044761)  
FIBER_1 = (24.446876,118.052855)

#依照时间筛选船只的AIS数据
# MMSI = 413788251
# Time_start = '19/07/24 09:48'
# Time_end = '19/07/24 09:53'
MMSI =  412440690
Time_start = '19/07/24 09:36:15'
Time_end = '19/07/24 09:36:30'


filepath = '/home/huangwj/DAS/BoatTrajectory/'
flag=1
if flag ==1:
    PosFile='AIS_DATA_RAW/gulangyu/pos_bj0829_20240715120000_20240816120000_3501.csv'
    StaticFile='AIS_DATA_RAW/gulangyu/static_bj0829_20240715120000_20240816120000_3501.csv'
    ST_UTC8=datetime.datetime.strptime("15/07/24 00:01", "%d/%m/%y %H:%M")
    ET_UTC8=datetime.datetime.strptime("15/08/24 23:59", "%d/%m/%y %H:%M")
    filename='AIS_SHIP_DATA/FiberBoatMessage_gly_240715_240815.csv'

dfp=pd.read_csv(PosFile,encoding='GBK',error_bad_lines =True )
df=pd.read_csv(StaticFile,encoding='GBK')
print(f'statis:{len(df)}',f'pos:{len(dfp)}')
ST_UTC0=ST_UTC8-timedelta(hours=8)
ET_UTC0=ET_UTC8-timedelta(hours=8)


dfp.rename(columns={'MMSI':'MMSI','经度':'lng','纬度':'lat','航向':'heading','航速(节)':'speed','航艏向':'Ship heading','时间':'Time'},inplace=True)

dfp = dfp[dfp['Time'].notnull()]
print(dfp.dtypes)
dfp['Time'] = pd.to_datetime(dfp['Time'],errors='ignore',format = '%Y-%m-%d %H:%M:%S')
dfp_0 = dfp[dfp['MMSI'] == MMSI]
print(f'statis:{len(df)}',f'pos:{len(dfp)}')

ts =datetime.datetime.strptime(Time_start, "%d/%m/%y %H:%M:%S")
te =datetime.datetime.strptime(Time_end, "%d/%m/%y %H:%M:%S")
dfp = dfp[dfp['Time']>=ts]
dfp = dfp[dfp['Time']<=te]
#%%
Ship_set = set(dfp['MMSI'].values)
Ship_set = [MMSI]
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
print(Ship_cross_fiber)
m = folium.Map(location=[24.4,118.04],width=1200,height=600,zoom_start=12)
folium.PolyLine([[FIBER_0[0],FIBER_0[1]],[FIBER_1[0],FIBER_1[1]]],color='green').add_to(m)
for ship in Ship_cross_fiber:
    dfp_ship = dfp[dfp['MMSI']==ship]
    dfp_ship.sort_values(by=['Time'],ascending = True, inplace = True)
    print(dfp_ship)
    lng_list = dfp_ship['lng'].tolist()
    lat_list = dfp_ship['lat'].tolist()
    for i in range(0,len(lng_list)-1):
        folium.PolyLine([[lat_list[i],lng_list[i]],[lat_list[i+1],lng_list[i+1]]]).add_to(m)
        folium.Circle([lat_list[i],lng_list[i]],color = 'red').add_to(m)
folium.Circle([lat_list[0],lng_list[0]],color = 'blue').add_to(m)
folium.Circle([lat_list[-1],lng_list[-1]],color = 'green').add_to(m)

print(f'AIS记录数量{len(dfp_ship)}')

display(m)
m.save('map.html')


#%%

# #读取csv
# os.chdir('/home/huangwj/DAS/DataSet/')
# POS_TIME = pd.read_csv('merge_pos_time.csv',index_col=0)
# POS_TIME['t0'] = pd.to_datetime(POS_TIME['t0'],errors='ignore',format = '%Y-%m-%d %H:%M:%S')
# POS_TIME['t1'] = pd.to_datetime(POS_TIME['t1'],errors='ignore',format = '%Y-%m-%d %H:%M:%S')
# POS_TIME['Delta_AIS_TIME'] = POS_TIME['t1']-POS_TIME['t0']

# #筛选时间间距较少的记录
# POS_TIME['Delta_AIS_TIME'] = POS_TIME['Delta_AIS_TIME'].apply(lambda x: x.total_seconds())
# POS_TIME = POS_TIME[POS_TIME['Delta_AIS_TIME']<=10]
# POS_cal = POS_TIME['POS'][POS_TIME['Speed'].astype(float)>8].tolist()

# #筛选不正常的数据点

# POS_cal = [pos for pos in POS_cal if eval(pos)[1]<22.169]





