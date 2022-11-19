
# %%
import numpy as np
import pandas as pd
import math
import datetime
from datetime import timedelta
from math import sin,cos



def direct(V1:np.array,V2:np.array):
    if (V1[0]*V2[1]-V1[1]*V2[0])<0:
        sign=-1
    else:
        sign=1
    return sign
def cross(A,B,C,D):  #判断线段AB 与CD是否相交
    cross=0  #不相交
    if direct(C-A,D-A)*direct(C-B,D-B)<0:
        if direct(A-C,B-C)*direct(A-D,B-D)<0:
            cross=1
    return cross
def cross_point(A,B,C,D):   ##  计算两直线的交点
    x1=A[0]
    y1=A[1]
    x2=B[0]
    y2=B[1]
    x3=C[0]
    y3=C[1]
    x4=D[0]
    y4=D[1]
    k1=(y2-y1)*1.0/(x2-x1)#计算k1,由于点均为整数，需要进行浮点数转化
    b1=y1*1.0-x1*k1*1.0#整型转浮点型是关键
    if (x4-x3)==0:#L2直线斜率不存在操作
        k2=None
        b2=0
    else:
        k2=(y4-y3)*1.0/(x4-x3)#斜率存在操作
        b2=y3*1.0-x3*k2*1.0
    if k2==None:
        x=x3
    else:
        x=(b2-b1)*1.0/(k1-k2)
    y=k1*x*1.0+b1*1.0
    #判断CD在地图上的朝向
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

    return [x,y],Tra_dirc

def GetCrossAngle(A,B,C,D):
    x1=A[0]
    y1=A[1]
    x2=B[0]
    y2=B[1]
    x3=C[0]
    y3=C[1]
    x4=D[0]
    y4=D[1]
    arr_0 = np.array([(x2 - x1), (y2 - y1)])
    arr_1 = np.array([(x4 - x3), (y4 - y3)])
    cos_value = (float(arr_0.dot(arr_1)) / (np.sqrt(arr_0.dot(arr_0)) * np.sqrt(arr_1.dot(arr_1))))
    if cos_value>1:
        cos_value=1
    elif cos_value<-1:
        cos_value=-1
    return np.arccos(cos_value)


EARTH_REDIUS = 6378.137
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

def crosstime(tra_0,tra_1,cp,time_0,time_1):
    s1=getDistance(cp[0],cp[1],tra_0[0],tra_0[1])
    S=getDistance(tra_1[0],tra_1[1],tra_0[0],tra_0[1])
    deltaT=(time_1-time_0)*s1/S
    t=time_0+deltaT
    return t

def crossSpeed(tra_0,tra_1,cp,speed_0,speed_1):
    s1=getDistance(cp[0],cp[1],tra_0[0],tra_0[1])
    S=getDistance(tra_1[0],tra_1[1],tra_0[0],tra_0[1])
    deltaS=(speed_1-speed_0)*s1/S
    t=speed_0+deltaS
    return t





def AISData(PosFile:str,StaticFile:str,ST_UTC8:datetime.datetime,ET_UTC8:datetime.datetime):
    df = pd.read_csv(PosFile) #UTC+8
    static=pd.read_csv(StaticFile)


    ST_UTC0=ST_UTC8-timedelta(hours=8)
    ET_UTC0=ET_UTC8-timedelta(hours=8)


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

    
    crossFiberBoat=pd.DataFrame(columns=['MMSI','CrossTime','Time_0','Time_1','Time_0_1(min)','CrossSpeed','Speed_0_1','lat','lng','disFromEnd/Km','tra_direction','Angle'])
    #Time_0_1(min)反映过光纤前船只坐标时间的差异，用以筛选掉过光纤前    后时间差异特别大的数据，默认5分钟内的数据比较合适。
    #Speed_0_1反映船只的速度变化，用来判断船只的状态是否不正常
    
    
    #遍历每一行
    df_time['MMSI'] = df_time['MMSI'].astype(str)
    INDEX=df_time['MMSI'].value_counts().index
    FIBER_0=np.array([22.140,113.709])
    FIBER_1=np.array([22.168,113.801])
    
    for mmsi in INDEX:
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
                [lati,lng],Tra_dirc=cross_point(FIBER_0,FIBER_1,    Tra_0,Tra_1)
                t=crosstime(Tra_0,Tra_1,[lati,lng],traTime[i],  traTime[i+1])
                s=crossSpeed(Tra_0,Tra_1,[lati,lng],speed[i],   speed[i+1])
                d=getDistance(FIBER_1[0], FIBER_1[1], lati, lng)
                crossFiberBoat.loc[len(crossFiberBoat.index)] =     (mmsi,t,traTime[i],traTime[i+1],(traTime[i+1]-traTime[i]).seconds/60,s, speed[i+1]-speed[i],lati,lng,d,Tra_dirc,math.    degrees(GetCrossAngle(FIBER_0,FIBER_1,Tra_0,    Tra_1)))
                
    crossFiberBoat['MMSI']=crossFiberBoat['MMSI'].astype('str')
    #删除过光纤前后时间差异较大的点，暂定为5分钟
    crossFiberBoat=crossFiberBoat[crossFiberBoat['Time_0_1(min)']<5]
    #将速度从节换算为m/s
    crossFiberBoat['CrossSpeed']=crossFiberBoat['CrossSpeed']*1000/3600*1.852
 
    crossFiberBoat.sort_values(by='CrossTime').to_excel ("crossFiberBoat.xlsx")
  
    #关联船的信息
    static['MMSI']=static['MMSI'].astype('str')
    MMSI=list(crossFiberBoat["MMSI"])
    MMSI=list(set(MMSI))
    
    FiberBoatMessage=pd.DataFrame(columns=['MMSI','CrossTime',  'Time_0','Time_1','CrossSpeed(m/s)','Delta_Speed','lat','lng','disFromEnd/Km', 'tra_direction','Angle','length','width','depth','type'])
    
    for mmsi in MMSI:
        tmp1=crossFiberBoat[crossFiberBoat['MMSI']==mmsi]
        tmp2=static[static['MMSI']==mmsi]
        tmp=np.array([tmp1['MMSI'],tmp1['CrossTime'],tmp1['Time_0'],tmp1['Time_1'],tmp1['CrossSpeed'],tmp1['Speed_0_1'],tmp1['lat'],tmp1['lng'],tmp1['disFromEnd/Km'],tmp1['tra_direction'],tmp1['Angle']])
        tmp=np.row_stack((tmp,tmp.shape[1]*list(tmp2["船长（米）"])))
        tmp=np.row_stack((tmp,tmp.shape[1]*list(tmp2['船宽（米）'])))
        tmp=np.row_stack((tmp,tmp.shape[1]*list(tmp2['吃水（米）'])))
        tmp=np.row_stack((tmp,tmp.shape[1]*list(tmp2['类型'])))
        df_temp=pd.DataFrame(np.transpose(tmp),columns=['MMSI', 'CrossTime',  'Time_0','Time_1','CrossSpeed(m/s)','Delta_Speed','lat','lng','disFromEnd/Km','tra_direction','Angle','length','width','depth',  'type']) #'CrossTime':datetime.datetime
    
        FiberBoatMessage = pd.concat([FiberBoatMessage,df_temp], ignore_index=False)
    
    
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

def AnchorShip(FiberBoatMessage,MINCHANNEL,n_channels,channel_spacing,ST,ET,TT):
    
    CT=list(FiberBoatMessage['CrossTime'])
    deltaCT=[]
    for i in range(0,len(CT)):
        deltaCT.append(round((CT[i]-ST)/(ET-ST)*TT))

    vline_indx = deltaCT

    #利用AIS数据在DAS数据上标定船只通过的时间以及位置
    Dist=list(FiberBoatMessage['disFromEnd/Km'])
    #区域距离标定(3.33和120都是经验参数)
    #RegionDistUpper=list(n_channels-(Dist+3.33)*1000/channel_spacing+120)
    RegionDistUpper=[round(n_channels-(i-MINCHANNEL+3.33)*1000/channel_spacing+200) for i in Dist]
    RegionDistdown=[round(n_channels-(i-MINCHANNEL+3.33)*1000/channel_spacing-200) for i in Dist]
    deltaCTUpper=[i+round(TT*0.05) for i in vline_indx]
    deltaCTdown=[i-round(TT*0.05) for i in vline_indx]
    CT=list(FiberBoatMessage['Time_0'])
    deltaCT=[]
    for i in range(0,len(CT)):
        deltaCT.append(round((CT[i]-ST)/(ET-ST)*TT))
    deltaCTdown=deltaCT
    CT=list(FiberBoatMessage['Time_1'])
    deltaCT=[]
    for i in range(0,len(CT)):
        deltaCT.append(round((CT[i]-ST)/(ET-ST)*TT))
    deltaCTUpper=deltaCT
    return deltaCTUpper,deltaCTdown,RegionDistUpper,RegionDistdown

def ShipTraj(PosFile:str,StaticFile:str,ST_UTC8:datetime.datetime,ET_UTC8:datetime.datetime,MMSI):#船只轨迹数据读取
    df = pd.read_csv(PosFile) #UTC+8
    static=pd.read_csv(StaticFile)
    ST_UTC0=ST_UTC8-timedelta(hours=8)
    ET_UTC0=ET_UTC8-timedelta(hours=8)
    df['MMSI'].value_counts()
    df['时间']=pd.to_datetime(df['时间'])
    df_time=pd.DataFrame(df[(df['时间']>ST_UTC8) &  (df['时间']<ET_UTC8)])
    return df_time[df_time['MMSI']==MMSI]

