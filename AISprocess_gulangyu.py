
# %%
import numpy as np
import pandas as pd
import math
import datetime
from datetime import timedelta
from math import sin,cos
import os
from tqdm import tqdm
from shapely.geometry import LineString
os.chdir('/home/huangwj/DAS/BoatTrajectory/')

EARTH_REDIUS = 6378.137

def cross(A,B,C,D):
    track1=[A,B]
    track2=[C,D]
    line1 = LineString(track1)
    line2 = LineString(track2)
    
    if line1.intersects(line2):
        return True
    else:
        return False


def cross_point(A,B,C,D):
    track1=[A,B]
    track2=[C,D]
    line1 = LineString(track1)    #经纬度坐标系
    line2 = LineString(track2)
    
    intersection = line1.intersection(line2)

        #判断CD在地图上的朝向
    x1=A[0]
    y1=A[1]
    x2=B[0]
    y2=B[1]
    x3=C[0]
    y3=C[1]
    x4=D[0]
    y4=D[1]
    LATD=['N','S']
    LNGD=['E','W']
    if (x4-x3)>0:
        direc_lat=[1,0]
    else: 
        direc_lat=[0,1]

    if (y4-y3)>0:
        direc_lng=[1,0]
    else: 
        direc_lng=[0,1]

    Tra_dirc=LNGD[direc_lng[0]]+LATD[direc_lat[0]]+'->'+LNGD[direc_lng[1]]+LATD[direc_lat[1]]
    
    if intersection.geom_type == 'Point':
        return [intersection.x, intersection.y],Tra_dirc
    elif intersection.geom_type == 'MultiPoint':
        # 如果有多个交点，可以选择其中一个或根据需求处理
        return [intersection[0].x, intersection[0].y],Tra_dirc
    else:
        return None



def convert_to_cartesian(latitude, longitude, radius):
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

def GetCrossAngle(A, B, C, D, radius=EARTH_REDIUS):
    """
    计算经纬度表示的两个线段之间的夹角
    
    参数：
    A, B, C, D: 包含经纬度信息的四个点坐标
    radius: 地球半径，单位为千米，默认值为6371千米
    
    返回值：
    两个线段之间的夹角，单位为弧度
    
    """
    # 将经纬度转换为直角坐标系
    point_A = convert_to_cartesian(A[0], A[1], radius)
    point_B = convert_to_cartesian(B[0], B[1], radius)
    point_C = convert_to_cartesian(C[0], C[1], radius)
    point_D = convert_to_cartesian(D[0], D[1], radius)
    
    # 计算两个向量之间的夹角
    vector_AB = point_B - point_A
    vector_CD = point_D - point_C
    
    cos_value = np.dot(vector_AB, vector_CD) / (np.linalg.norm(vector_AB) * np.linalg.norm(vector_CD))
    if cos_value > 1:
        cos_value = 1
    elif cos_value < -1:
        cos_value = -1
    
    return np.arccos(cos_value)



def rad(d):
    return d * math.pi / 180.0

def getDistance(lat1, lng1, lat2, lng2):
    radLat1 = rad(lat1)
    radLat2 = rad(lat2)
    a = radLat1 - radLat2
    b = rad(lng1) - rad(lng2)
    s = 2 * math.asin(math.sqrt(math.pow(sin(a/2), 2) + cos(radLat1) * cos(radLat2) * math.pow(sin(b/2), 2)))
    s = s * EARTH_REDIUS
    return s

def crosstime(tra_0,tra_1,cp,time_0,time_1,speed_0,speed_1):
    l=getDistance(cp[0],cp[1],tra_0[0],tra_0[1])
    L=getDistance(tra_1[0],tra_1[1],tra_0[0],tra_0[1])


    T=2*L/(speed_0+speed_1) #总耗时
    a=(speed_1-speed_0)/T #加速度
    v=math.sqrt(2*a*l+speed_0*speed_0)
    if a==0:
        deltaT=(time_1-time_0)*l/L
    else:
        deltaT=(time_1-time_0)*((v-speed_0)/a)/T

    """
    if pd.isna(deltaT):
        print(deltaT)
        print(time_1,time_0,speed_0,speed_1,l,L,T,a,v,deltaT)"""
    
    if np.isnan(a):   #当加速度值为nan时候，将船只移动过程认为是匀速运动 
        deltaT=(time_1-time_0)*l/L

    seconds = int(deltaT.total_seconds())
    deltaT = timedelta(seconds=seconds)
    t=time_0+deltaT

    return t

def crossSpeed(tra_0,tra_1,cp,speed_0,speed_1):
    l=getDistance(cp[0],cp[1],tra_0[0],tra_0[1])
    L=getDistance(tra_1[0],tra_1[1],tra_0[0],tra_0[1])
    '''
    deltaS=(speed_1-speed_0)*s1/S
    v=speed_0+deltaS
    '''
    T=2*L/(speed_0+speed_1) #总耗时
    a=(speed_1-speed_0)/T #加速度
    v=math.sqrt(2*a*l+speed_0*speed_0)
    return v

def average_speed(row):
    if row['Time_0_1(min)']>0:
        speed = getDistance(row['lat1'], row['lng1'], row['lat_0'], row['lng_0'])/row['Time_0_1(min)']/60
    else:
        speed = row['CrossSpeed']
    return speed*1000



def AISData(Pos_pd:pd.DataFrame,Static_pd:pd.DataFrame,ST_UTC8:datetime.datetime,ET_UTC8:datetime.datetime,FIBER_0=(22.161561, 113.712681), FIBER_1=(22.167286, 113.795594)):

    Pos_pd['MMSI'].value_counts()
    Pos_pd['时间']=pd.to_datetime(Pos_pd['时间'])
    Pos_pd=pd.DataFrame(Pos_pd[(Pos_pd['时间']>ST_UTC8) &  (Pos_pd['时间']<ET_UTC8)])

    INDEX=Pos_pd['MMSI'].value_counts().index
    temp=Pos_pd['MMSI'].value_counts()

    #删除只有一个数据的轨迹
    dropindex=[]
    for i in range(0,len(temp)):
        if temp.loc[INDEX[i]]==1:
            dropindex.append(temp.index[i])
    Pos_pd=Pos_pd[~Pos_pd['MMSI'].isin(dropindex)]

    #删除船速为0的数据
    Pos_pd=Pos_pd[Pos_pd['航速(节)']!=0]
    #print('Data clean!')    
    
    crossFiberBoat=pd.DataFrame(columns=['MMSI','CrossTime','Time_0','Time_1','Time_0_1(min)','CrossSpeed','Speed_0_1','lat','lng','lat_0','lng_0','lat1','lng1','disFromEnd/Km','tra_direction','Angle'])
    #Time_0_1(min)反映过光纤前船只坐标时间的差异，用以筛选掉过光纤前后时间差异特别大的数据，默认5分钟内的数据比较合适。
    #Speed_0_1反映船只的速度变化，用来判断船只的状态是否不正常
    
    #遍历每一行
    Pos_pd.loc[:,'MMSI'] = Pos_pd['MMSI'].astype(str)
    INDEX=Pos_pd['MMSI'].value_counts().index

    '''如果仅考虑直线段光纤，应该分析fiber_0到Max_distance通道间的das数据，约八公里长度的光纤数据'''
    print('计算过光纤船只的速度与方向')
    for mmsi in tqdm(INDEX):
        df=Pos_pd[Pos_pd['MMSI']==mmsi]
        long=df['经度']
        lat=df['纬度']
        traTime=df['时间']
        speed=df['航速(节)']
        long=list(long)
        lat=list(lat)
        traTime=list(traTime)
        speed=list(speed)
        #lng越大越往东，lat越大越往北
        for i in range(0,len(long)-1):
            Tra_0=np.array([lat[i],long[i]])
            Tra_1=np.array([lat[i+1],long[i+1]])
            if cross(FIBER_0,FIBER_1,Tra_0,Tra_1)==1:
                [lati,lng],Tra_dirc=cross_point(FIBER_0,FIBER_1,Tra_0,Tra_1)
                t=crosstime(Tra_0,Tra_1,[lati,lng],traTime[i], traTime[i+1],speed[i], speed[i+1])
                s=crossSpeed(Tra_0,Tra_1,[lati,lng],speed[i], speed[i+1])
                d=getDistance(FIBER_1[0], FIBER_1[1], lati, lng)
                crossFiberBoat.loc[len(crossFiberBoat.index)] = (mmsi,t,traTime[i],traTime[i+1],(traTime[i+1]-traTime[i]).total_seconds()/60,s, speed[i+1]-speed[i],lati,lng,Tra_0[0],Tra_0[1],Tra_1[0],Tra_1[1],d,Tra_dirc,math.degrees(GetCrossAngle(FIBER_0,FIBER_1,Tra_0,Tra_1)))
    crossFiberBoat.loc[:,'MMSI']=crossFiberBoat['MMSI'].astype('str')

    #将速度从节换算为m/s
    crossFiberBoat.loc[:,'CrossSpeed']=crossFiberBoat['CrossSpeed']*1000/3600*1.852
    #关联船的信息
    Static_pd.loc[:,'MMSI']=Static_pd['MMSI'].astype('str')
    MMSI=list(crossFiberBoat["MMSI"])
    MMSI=list(set(MMSI))
    FiberBoatMessage= pd.merge(crossFiberBoat, Static_pd, left_on='MMSI', right_on='MMSI',how='left')
    FiberBoatMessage.loc[:,'cross_position'] = FiberBoatMessage['disFromEnd/Km']  #根据光线重定位项目
    FiberBoatMessage.loc[:,"CrossTime"]=pd.to_datetime(FiberBoatMessage["CrossTime"])
    FiberBoatMessage.loc[:,"Time_0"]=pd.to_datetime(FiberBoatMessage["Time_0"])
    FiberBoatMessage.loc[:,"Time_1"]=pd.to_datetime(FiberBoatMessage["Time_1"])
    FiberBoatMessage['average_speed'] = FiberBoatMessage.apply(lambda row:average_speed(row),axis=1) #计算平均速度
    FiberBoatMessage['err_speed'] = (FiberBoatMessage['average_speed']-FiberBoatMessage['CrossSpeed']).abs()/FiberBoatMessage['CrossSpeed']
    FiberBoatMessage = FiberBoatMessage[FiberBoatMessage['Time_0_1(min)']<10]
    return FiberBoatMessage
# %%
if __name__=='__main__':
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
    print(df.columns)
    print(dfp.columns)

    ST_UTC0=ST_UTC8-timedelta(hours=8)
    ET_UTC0=ET_UTC8-timedelta(hours=8)

    #鼓浪屿
    # FIBER_0 = (24.454861,118.041744)
    # FIBER_1 = (24.447322,118.053203)
    FIBER_0 = (24.454581,118.044761)  
    FIBER_1 = (24.446876,118.052855)
    
    
    
    FIBER_0 = (24.454556,118.044845)  #重定位后经纬度
    FIBER_1 = (24.446949,118.052813)

    print(len(dfp),len(df))

    FiberBoatMessage=AISData(dfp,df,ST_UTC8,ET_UTC8,FIBER_0,FIBER_1)
    FiberBoatMessage.sort_values(by='CrossTime',ascending=True,inplace=True,ignore_index=True)
    print(len(FiberBoatMessage))
    #修改FiberBoatMessage的列名
    FiberBoatMessage.rename(columns = {'船名' : 'ShipName', '船旗' : 'Flag','船长（米）':'length','船宽（米）':'width','吃水（米）':'depth','类型':'type','更新时间':'UpdateTime'}, inplace = True)
    

    FiberBoatMessage.to_csv(filename)

    print(f'Write {filename}')
