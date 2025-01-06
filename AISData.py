
# %%
import numpy as np
import pandas as pd
import math
import datetime
from datetime import timedelta
from math import sin,cos,radians,acos,sqrt,degrees,atan2
from tqdm import tqdm
from shapely.geometry import LineString


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





def AISData(PosDf:pd.DataFrame,StaticDF:pd.DataFrame,ST_UTC8:datetime.datetime,ET_UTC8:datetime.datetime,FIBER_0=(22.161561, 113.712681), FIBER_1=(22.167286, 113.795594)):
    #FIBER_0和FIBER_1默认为三角岛和桂山岛之间光缆起止点的经纬度
    df = PosDf #UTC+8
    static=StaticDF
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

    #删除船速为0的数据
    df_time=df_time[df_time['航速(节)']!=0]
    print('Data clean!')    
    crossFiberBoat=pd.DataFrame(columns=['MMSI','CrossTime','Time_0','Time_1','Time_0_1(min)','CrossSpeed','Speed_0_1','lat','lng','disFromEnd/Km','tra_direction','Angle'])
    #Time_0_1(min)反映过光纤前船只坐标时间的差异，用以筛选掉过光纤前    后时间差异特别大的数据，默认5分钟内的数据比较合适。
    #Speed_0_1反映船只的速度变化，用来判断船只的状态是否不正常
    
    #遍历每一行
    df_time['MMSI'] = df_time['MMSI'].astype(str)
    INDEX=df_time['MMSI'].value_counts().index
    # FIBER_0=(22.161561, 113.712681)  #纬度（latitude）  经度  （longtitude） （调整后光纤位置）
    # FIBER_1=(22.167286, 113.795594) #桂山岛方向 

    '''如果仅考虑直线段光纤，应该分析fiber_0到Max_distance通道间的das数据，约八公里长度的光纤数据'''
    print('计算过光纤船只的速度与方向')
    for mmsi in tqdm(INDEX):
        df=df_time[df_time['MMSI']==mmsi]
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
                crossFiberBoat.loc[len(crossFiberBoat.index)] = (mmsi,t,traTime[i],traTime[i+1],(traTime[i+1]-traTime[i]).total_seconds()/60,s, speed[i+1]-speed[i],lati,lng,d,Tra_dirc,math.degrees(GetCrossAngle(FIBER_0,FIBER_1,Tra_0,Tra_1)))
                
    crossFiberBoat['MMSI']=crossFiberBoat['MMSI'].astype('str')
    #删除过光纤前后时间差异较大的点，暂定为6分钟
    crossFiberBoat=crossFiberBoat[crossFiberBoat['Time_0_1(min)']<6]
    #将速度从节换算为m/s
    crossFiberBoat['CrossSpeed']=crossFiberBoat['CrossSpeed']*1000/3600*1.852
 
    #crossFiberBoat.sort_values(by='CrossTime').to_excel ("crossFiberBoat.xlsx")
  
    #关联船的信息
    static['MMSI']=static['MMSI'].astype('str')
    MMSI=list(crossFiberBoat["MMSI"])
    MMSI=list(set(MMSI))
    
    #FiberBoatMessage=pd.DataFrame(columns=['MMSI','CrossTime',  'Time_0','Time_1','CrossSpeed(m/s)','Delta_Speed','lat','lng','disFromEnd/Km', 'tra_direction','Angle','length','width','depth','type'])
    print('关联船只信息')

    FiberBoatMessage= pd.merge(crossFiberBoat, static, left_on='MMSI', right_on='MMSI',how='left')
    FiberBoatMessage["CrossTime"]=pd.to_datetime(FiberBoatMessage["CrossTime"])
    FiberBoatMessage["Time_0"]=pd.to_datetime(FiberBoatMessage["Time_0"])
    FiberBoatMessage["Time_1"]=pd.to_datetime(FiberBoatMessage["Time_1"])
    return FiberBoatMessage
# %%
#返回筛选结果，以及每个筛选条件下船只的数量，以供调整筛选条件
def FilterBoat(FiberBoatMessage:pd.DataFrame,DT:datetime.datetime,Angle:list,Dist:list,Speed:list):
    FiberBoatMessage['CrossTime'].iloc[1]
    index=[]
    for i in range(1,len(FiberBoatMessage['CrossTime'])-1):
        t0=FiberBoatMessage['CrossTime'].iloc[i-1]
        t1=FiberBoatMessage['CrossTime'].iloc[i]
        t2=FiberBoatMessage['CrossTime'].iloc[i+1]
        if (t1-t0>=DT) and (t2-t1>=DT):
            index.append(i)
    Boat_TimeLimit=FiberBoatMessage.iloc[index]

    #筛选垂直过光纤的船只
    Boat_angel=Boat_TimeLimit[Boat_TimeLimit['Angle']>Angle[0]]
    Boat_angel=Boat_angel[Boat_angel['Angle']<Angle[1]] 

    #筛选船只经过海中心
    Boat_dis=Boat_angel[Boat_angel['disFromEnd/Km']>Dist[0]] 
    Boat_dis=Boat_dis[Boat_dis['disFromEnd/Km']<Dist[1]] 
    #筛选高速船
    Boat_speed=Boat_dis[Boat_dis['CrossSpeed(m/s)']>Speed[0]]
    Boat_speed=Boat_speed[Boat_speed['CrossSpeed(m/s)']<Speed[1]]
    return Boat_speed,[len(Boat_TimeLimit),len(Boat_angel),len(Boat_dis),len(Boat_speed)]

def cal_channel_by(depth, zero_offset, channel_spacing):
    # 根据depth反推最接近的channel
    n_channels_ceil = np.ceil((depth - zero_offset) / channel_spacing) + 1
    n_channels_floor = np.floor((depth - zero_offset) / channel_spacing) + 1
    if (depth <= zero_offset + (n_channels_ceil - 1) * channel_spacing):
        return int(n_channels_ceil)
    else:
        return int(n_channels_floor)
'''
def AnchorShip(FiberBoatMessage,MINCHANNEL,MAXCHANNEL,n_channels,channel_spacing,zero_offset,ST,ET,TT):
    CT=list(FiberBoatMessage['CrossTime'])
    deltaCT=[]
    for i in range(0,len(CT)):
        deltaCT.append(round((CT[i]-ST).total_seconds()/(ET-ST).total_seconds()*TT))
    vline_indx = deltaCT
    Fiber_1 = 13.07 #Km 桂山岛起点，用于确定cross_pos
    #利用AIS数据在DAS数据上标定船只通过的时间以及位置
    #Dist=list(FiberBoatMessage['disFromEnd/Km'])
    Dist = list(Fiber_1 - np.array(FiberBoatMessage['disFromEnd/Km']))  #根据光线重定位项目

    #区域距离标定(3.33和120都是经验参数)
    #RegionDistUpper=list(n_channels-(Dist+3.33)*1000/channel_spacing+120)
    #RegionDistUpper=[round(n_channels-(i-MINCHANNEL+3.33)*1000/channel_spacing+200) for i in Dist]
    #RegionDistdown=[round(n_channels-(i-MINCHANNEL+3.33)*1000/channel_spacing-200) for i in Dist]
    #ChannelOffSet=round(3.33*1000/channel_spacing)
    #DISTFROMEND=[n_channels-round(i*1000/channel_spacing)-ChannelOffSet for i in Dist]

    DISTFROMEND = []
    for dist in Dist:
        DISTFROMEND.append(cal_channel_by(dist*1000, zero_offset, channel_spacing))
    MINCHANNEL = cal_channel_by(MINCHANNEL*1000, zero_offset, channel_spacing)
    MAXCHANNEL = cal_channel_by(MAXCHANNEL*1000, zero_offset, channel_spacing)

    RegionDistUpper=[]
    RegionDistdown=[]
    RegionArea=200 #channel
    ShipIndex=[]
    for dist in DISTFROMEND:
        if dist>MINCHANNEL and dist<MAXCHANNEL:
            RegionDistUpper.append(min(dist+RegionArea,MAXCHANNEL))
            RegionDistdown.append(max(dist-RegionArea,MINCHANNEL))
            ShipIndex.append(DISTFROMEND.index(dist))
    RegionDistUpper=list(np.array(RegionDistUpper)-MINCHANNEL)   
    RegionDistdown=list(np.array(RegionDistdown)-MINCHANNEL)   
    

    deltaCTUpper=[i+round(TT*0.05) for i in vline_indx]
    deltaCTdown=[i-round(TT*0.05) for i in vline_indx]
    CT=list(FiberBoatMessage['Time_0'])
    deltaCT=[]
    for i in range(0,len(CT)):
        deltaCT.append(round((CT[i]-ST)/(ET-ST)*TT))
    deltaCTdown=[]
    for i in ShipIndex:
        deltaCTdown.append(deltaCT[i])


    CT=list(FiberBoatMessage['Time_1'])
    deltaCT=[]
    for i in range(0,len(CT)):
        deltaCT.append(round((CT[i]-ST)/(ET-ST)*TT))
    deltaCTUpper=[]   
    for i in ShipIndex:
        deltaCTUpper.append(deltaCT[i])
    
    return deltaCTUpper,deltaCTdown,RegionDistUpper,RegionDistdown
'''
def cal_time(t0,t1,ST,ET):
    T0 = t0
    T1 = t1
    if t0<ST:
        T0 = ST
    if t1>ET:
        T1 = ET
    return T0,T1



def AnchorShip(FiberBoatMessage,MINCHANNEL,MAXCHANNEL,n_channels,channel_spacing,zero_offset,ST,ET,Time_dimension,Dist_dimension):
    deltaCTUpper = []
    deltaCTdown = []
    RegionDistUpper = []
    RegionDistdown = []
    MMSI = []
    TotalSeconds = (ET-ST).total_seconds()
    minchannel = cal_channel_by(MINCHANNEL*1000, zero_offset, channel_spacing)
    maxchannel = cal_channel_by(MAXCHANNEL*1000, zero_offset, channel_spacing)
    RegionArea=200 #channel
    for row_index in range(0,len(FiberBoatMessage)):
        row = FiberBoatMessage.iloc[row_index,:]
        Fiber_1 = 13.07 #Km 桂山岛起点，用于确定cross_pos
        #利用时间和位置约束来筛选船舶数据
        if (row['CrossTime']>=ST) and (row['CrossTime']<=ET):
            Dist = Fiber_1 - row['disFromEnd/Km'] #根据光线重定位项目
            if (Dist>= MINCHANNEL) and (Dist<= MAXCHANNEL):
                t0,t1 = cal_time(row['Time_0'],row['Time_1'],ST,ET)
                deltaCTdown.append(round((t0-ST).total_seconds()/TotalSeconds*Time_dimension))
                deltaCTUpper.append(round((t1-ST).total_seconds()/TotalSeconds*Time_dimension))
                dist = cal_channel_by(Dist*1000, zero_offset, channel_spacing)
                RegionDistUpper.append(min(dist+RegionArea,maxchannel)-minchannel)
                RegionDistdown.append(max(dist-RegionArea,minchannel)-minchannel)
                RegionDistUpper.append(round((min(dist+RegionArea,maxchannel)-minchannel)/(maxchannel-minchannel)*Dist_dimension))
                RegionDistdown.append(round((max(dist-RegionArea,minchannel)-minchannel)/(maxchannel-minchannel)*Dist_dimension))
                MMSI.append(row['MMSI'])
    return deltaCTUpper,deltaCTdown,RegionDistUpper,RegionDistdown,MMSI

def ShipTraj(PosFile:str,StaticFile:str,ST_UTC8:datetime.datetime,ET_UTC8:datetime.datetime,MMSI):#船只轨迹数据读取
    df = pd.read_csv(PosFile) #UTC+8
    static=pd.read_csv(StaticFile)
    ST_UTC0=ST_UTC8-timedelta(hours=8)
    ET_UTC0=ET_UTC8-timedelta(hours=8)
    df['MMSI'].value_counts()
    df['时间']=pd.to_datetime(df['时间'])
    df_time=pd.DataFrame(df[(df['时间']>ST_UTC8) &  (df['时间']<ET_UTC8)])
    return df_time[df_time['MMSI']==MMSI]


import glob
import re
from datetime import datetime
def Index_AIS_csv_by_DAS_Data(input_date:str):
    file_pattern = '/home/huangwj/DAS/BoatTrajectory/AIS_SHIP_DATA/FiberBoatMessage*.csv'
    # 获取匹配文件路径列表
    file_list = glob.glob(file_pattern)
    matched_files=[]
    # 正则表达式模式用于从文件名中提取日期
    date_pattern = r'FiberBoatMessage_(\d{6})_(\d{6})\.csv'
    for file_path in file_list:
        # 提取文件名中的日期信息
        match = re.search(date_pattern, file_path)
        if match:
            file_start_date = datetime.strptime(match.group(1), '%y%m%d')
            file_end_date = datetime.strptime(match.group(2), '%y%m%d')
             # 提取日期信息，添加到列表中
            if file_start_date <= input_date <= file_end_date:
                matched_files.append(file_path)
    return matched_files


if __name__=='__main__':
    # from DASFileRead import DasFileRead

    # SHIP=0
    # MMSI=['413260090','413208430','413471740','413226010','413231470','413260090']
    # DataPath='/home/huangwj/DAS/BoatTrajectory/DataforAIS'
    # StartTime=["24/07/22 09:54","24/07/22 10:01","24/07/22 21:09",'10/09/21 13:17','06/09/21 11:56','09/09/21 9:33']
    # EndTime=["24/07/22 10:00","24/07/22 10:03","24/07/22 21:11",'10/09/21 13:22','06/09/21 11:58','09/09/21 9:37']

    # MINTIME = [0.5,0,0,0,0,0]   #0.5  0
    # MAXTIME = [-1,-1,-1,-1,2.5,-1]   #2    3
    # MINCHANNEL=[7.8,10,10,11,10,9]   #8.5   7.8  
    # MAXCHANNEL=[10.5,13,13,13.5,12.5,10.5] #Km  #9.7  10.5
    # MINTIME=MINTIME[SHIP]
    # MAXTIME=MAXTIME[SHIP]
    # MINCHANNEL=MINCHANNEL[SHIP]
    # MAXCHANNEL=MAXCHANNEL[SHIP]

    # ST_UTC8=datetime.strptime(StartTime[SHIP], "%d/%m/%y %H:%M")
    # ET_UTC8=datetime.strptime(EndTime[SHIP], "%d/%m/%y %H:%M")
    # ST_UTC0=ST_UTC8-timedelta(hours=8)
    # ET_UTC0=ET_UTC8-timedelta(hours=8)
    # FileSet,times=DasFileRead(ST_UTC0,ET_UTC0,DataPath)

    # ST=times[0]
    # ET=times[-1]+timedelta(minutes=1)
    # ST=ST+timedelta(hours=8)
    # ET=ET+timedelta(hours=8)


    # AIS_file_list = Index_AIS_csv_by_DAS_Data(datetime.strptime(StartTime[SHIP], "%d/%m/%y %H:%M"))
    # print(AIS_file_list)
    # FiberBoatMessage = pd.read_csv(AIS_file_list[0])
    # FiberBoatMessage['CrossTime'] = pd.to_datetime(FiberBoatMessage['CrossTime'],format="%Y-%m-%d %H:%M:%S")
    # FiberBoatMessage['Time_1'] = pd.to_datetime(FiberBoatMessage['Time_1'],format="%Y-%m-%d %H:%M:%S")
    # FiberBoatMessage['Time_0'] = pd.to_datetime(FiberBoatMessage['Time_0'],format="%Y-%m-%d %H:%M:%S")

    # n_channels = 4096
    # channel_spacing = 4.0838
    # zero_offset = 1.1048660452254646
    # deltaCTUpper,deltaCTdown,RegionDistUpper,RegionDistdown,MMSI = AnchorShip(FiberBoatMessage,MINCHANNEL,MAXCHANNEL,n_channels,channel_spacing,zero_offset,ST,ET,1500,600)
    # print(deltaCTUpper)
    # print(RegionDistdown,RegionDistUpper)
    # print(MMSI)
    
    FIBER_0 = (24.454861,118.041744)
    FIBER_1 = (24.447322,118.053203)
    #print(getDistance(FIBER_0[0],FIBER_0[1],FIBER_1[0],FIBER_1[1]))
