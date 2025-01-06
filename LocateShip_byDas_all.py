from AISData  import *
import datetime
import folium
from IPython.display import display
import matplotlib.pyplot as plt
from numpy import NaN

def calCrossPoint(pos0,pos1,sp0,sp1,t0:datetime.datetime,t1:datetime.datetime,ct:datetime.datetime):
    #计算pos0-pos1的距离
    L=getDistance(pos0[0], pos0[1],pos1[0], pos1[1])*1000
    #计算船舶加速度
    acc = (sp1-sp0)/(t1-t0).total_seconds()
    #在通行前，船舶运动的距离
    l = ((ct-t0).total_seconds())*(sp0+(sp0+acc*(ct-t0).total_seconds()))/2
    #计算通行事件的经纬度信息
    lat = pos0[0]+(pos1[0]-pos0[0])*(l/L)
    lng = pos0[1]+(pos1[1]-pos0[1])*(l/L)
    return lat,lng


def getAIS_byCrossTimeDas(mmsi:str,Ct:datetime.datetime,ais_df:pd.DataFrame):
    dT = datetime.timedelta(minutes=5)
    mmsi_df = ais_df[ais_df['MMSI']==mmsi]
    mmsi_df = mmsi_df[(mmsi_df['时间']<=(Ct+dT))&(mmsi_df['时间']>=(Ct-dT))]
    mmsi_df.sort_values(by='时间',ascending=True,inplace=True)
    if (len(mmsi_df[(mmsi_df['时间']<=Ct)])>=1)&(len(mmsi_df[(mmsi_df['时间']>=Ct)])>=1):
        lat0,lng0,sp0,t0 = mmsi_df[(mmsi_df['时间']<=Ct)][['纬度','经度','航速(节)','时间']].iloc[-1] #通行前的AIS信息
        lat1,lng1,sp1,t1 = mmsi_df[(mmsi_df['时间']>=Ct)][['纬度','经度','航速(节)','时间']].iloc[0]  #通行后的AIS信息
        sp0 = sp0*1000/3600*1.852
        sp1 = sp1*1000/3600*1.852
        lat,lng = calCrossPoint((lat0,lng0),(lat1,lng1),sp0,sp1,t0,t1,Ct)
        return mmsi_df,lat,lng
    else:
        return NaN,NaN,NaN


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



def Reposition2D(POS_cal:pd.DataFrame,FIBER0,FIBER1):
    
    x = POS_cal['lat'].to_numpy()
    y = POS_cal['lng'].to_numpy()
    n = len(x)
    sum_x = np.sum(x)
    sum_y = np.sum(y)
    sum_xy = np.sum(x * y)
    sum_x2 = np.sum(x**2)
    a = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x**2)
    b = (sum_y - a * sum_x) / n
    
    line_lat = np.linspace(FIBER0[0],FIBER1[0],50)
    line_lng = a*line_lat+b
    return (line_lat[0],line_lng[0]),(line_lat[-1],line_lng[-1]),a,b





def Reposition(POS_cal,FIBER0,FIBER1):
    x_c = []
    y_c = []
    z_c = []
    for i in range(0,len(POS_cal)):
        record = POS_cal.iloc[i]
        lng= record['lng']
        lat= record['lat']
        x,y,z = convert_to_cartesian(lat, lng)
        x_c.append(x)
        y_c.append(y)
        z_c.append(z)
    # %%在三维笛卡尔坐标系下线性拟合直线
    samples = np.array([x_c,y_c,z_c])
    A=np.array([np.mean(x_c),np.mean(y_c),np.mean(z_c)]).reshape(3,1)
    Y = samples-A
    s= np.zeros((3,3))
    for i in range(0,Y.shape[1]):
        y_vec= Y[:,0].reshape(3,1)
        y_vec_t = y_vec.T
        yty= y_vec_t @ y_vec *np.eye(3)
        yyt= y_vec @ y_vec_t
        s = s+(yty-yyt)
    eig_vals, eig_vecs = np.linalg.eig(s)
    idx = eig_vals.argsort() 
    # 最小的特征值
    minval = eig_vals[idx[0]]
    # 最小特征值对应的特征向量
    minvec = eig_vecs[:,idx[0]]
    point = A
    # 直线的方向向量
    direction = np.reshape(minvec,(3,1))
    # 生成直线上的点
    t = np.linspace(-7.018, 1.47, 100).reshape(1,-1)
    line_points = point + t[:,np.newaxis] * direction
    line_points=line_points[0,:,:]
    # 直线拟合三维图绘制
    fig = plt.figure(dpi=300)
    ax = plt.axes(projection='3d')
    ax.plot(line_points[0,:], line_points[1,:],line_points[2,:], 'orange')  
    ax.scatter(point[0], point[1], point[2], c='red')
    ax.scatter(x_c, y_c, z_c, c='blue')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title('Position fitting of submarine optical cables \nin a three-dimensional Cartesian coordinate system',fontsize = 8)
    plt.savefig('estimated_line.png', bbox_inches='tight', pad_inches=0.1, dpi=300)
    plt.close()

    line_lng=[]
    line_lat=[]
    for i in range(0,line_points.shape[1]):
        vec = line_points[:,i]
        lng,lat = cartesian_to_geo(vec[0], vec[1], vec[2])
        line_lng.append(lng)
        line_lat.append(lat)

    m = folium.Map(location=[24.4,118.04],width=1200,height=600,zoom_start=12)
    for i in range(0,len(POS_cal)):
        record = POS_cal.iloc[i]
        lng= record['lng']
        lat= record['lat']
        folium.CircleMarker(location = [lat, lng],radius=0.5,color = 'red').add_to(m)

    folium.CircleMarker(location = [line_lat[0],line_lng[0]],radius=0.5,color = 'Yellow').add_to(m)
    folium.CircleMarker(location = [line_lat[-1],line_lng[-1]],radius=0.5,color = 'Blue').add_to(m)
    #folium.PolyLine([[FIBER_0[0],FIBER_0[1]],[FIBER_1[0],FIBER_1[1]]],color='black').add_to(m)
    folium.PolyLine([[FIBER0[0],FIBER0[1]],[FIBER1[0],FIBER1[1]]],weight =0.5,color='Green').add_to(m)
    folium.PolyLine([[line_lat[0],line_lng[0]],[line_lat[-1],line_lng[-1]]],weight =0.5,color='orange').add_to(m)
    m.save('map_revise.html')
    return (line_lat[0],line_lng[0]),(line_lat[-1],line_lng[-1])

def calculate_distance_and_intersection(m, b, x0, y0):


    # 计算点到直线的距离
    distance = abs(m * x0 - y0 + b) / (m**2 + 1)**0.5

    # 计算垂直于直线的交点
    # 垂线的斜率是原直线斜率的负倒数
    m_perpendicular = -1 / m
    # 垂线通过点P(x0, y0)，使用点斜式方程
    b_perpendicular = y0 - m_perpendicular * x0
    # 解直线方程组找到交点
    x_intersection = (b_perpendicular - b) / (m - m_perpendicular)
    y_intersection = m * x_intersection + b

    return distance, (x_intersection, y_intersection)


# def fiber_end(lat,lng,lat0,lng0,Das_km,Fiber_end_km):
#     return getDistance(lat0,lng0, lat,lng) - abs(Das_km-Fiber_end_km)

def fiber_end(point0,point,Das_km,Fiber_end_km):
    return getDistance(point0[0],point0[1],point[0],point[1]) - abs(Das_km-Fiber_end_km)


def bisection_method( lat0,Das_km,Fiber_end_km,cable_slope,cable_bias, tol=1e-6):
    # 检查f(a)和f(b)是否异号
    if Fiber_end_km>=Das_km: #lat decrease,lng increase
        delta_lat = -0.2
    else:
        delta_lat = 0.2
    POINT = (lat0,cable_slope*(lat0)+cable_bias)
    point0 = POINT
    point = (lat0+delta_lat,cable_slope*(lat0+delta_lat)+cable_bias)

    if fiber_end(POINT,point0,Das_km,Fiber_end_km) * fiber_end(POINT,point,Das_km,Fiber_end_km) > 0:
        raise ValueError("f(a)和f(b)必须异号")
    if Fiber_end_km>=Das_km: 
        while (point0[0] - point[0]) / 2 > tol:
            point1 = ((point[0] + point0[0]) / 2, cable_slope*((point[0] + point0[0]) / 2)+cable_bias)
            if fiber_end(POINT,point1,Das_km,Fiber_end_km) == 0:
                return point1
            elif fiber_end(POINT,point0,Das_km,Fiber_end_km) * fiber_end(POINT,point1,Das_km,Fiber_end_km) < 0:
                point = point1
            else:
                point0 = point1
    else:
        while (point[0] - point0[0]) / 2 > tol:
            point1 = ((point[0] + point0[0]) / 2, cable_slope*((point[0] + point0[0]) / 2)+cable_bias)
            if fiber_end(POINT,point1,Das_km,Fiber_end_km) == 0:
                return point1
            elif fiber_end(POINT,point0,Das_km,Fiber_end_km) * fiber_end(POINT,point1,Das_km,Fiber_end_km) < 0:
                point = point1
            else:
                point0 = point1
                
            
    print(point[0],point0[0])
    return (point[0] + point0[0]) / 2




if __name__ =='__main__':
    #在海图中手动拾取的经纬度信息
    #FIBER_4---FIBER_3---FIBER_2---FIBER_0---FIBER_1---FIBER_5
    #FIBER0--FIBER1 (2.7km-3.85km)
    FIBER_0 = (24.454581,118.044761)  
    FIBER_1 = (24.446876,118.052855)
    FIBER_2 = (24.454451,118.043477)  #施工图中坐标
    FIBER_3 = (24.455283,118.042261)
    FIBER_4 = (24.454975,118.041261)
    FIBER_5 = (24.447589,118.054596)
        
    print('海底光缆的长度(Fiber0-Fiber1):',getDistance(FIBER_0[0],FIBER_0[1],FIBER_1[0],FIBER_1[1],)*1000)
    #在DAS数据中手动标注的船舶通行事件
    SHIP = pd.read_csv('/home/huangwj/DAS/BoatTrajectory/locate_fiber_data/ship_pick.csv')

    # SHIP = pd.DataFrame(columns=['MMSI',"CrossTimeDas",'CrossDistDas','Direction','lat','lng']) #columns=['MMSI',"CrossTimeDas",'CrossDistDas','Direction','lat','lng']
    # SHIP['MMSI'] = ['413788251','412440943','413788251','413788251','412443431','413875287','412440943','412440943','412440690']
    # SHIP["CrossTimeDas"] = ['2024/07/19 09:50:17.5','2024/08/08 15:29:09','2024/08/14 00:46:29','2024/08/06 12:07:41.5','2024/08/10 15:44:41.5','2024/7/18 15:28:56','2024/8/5 09:33:09','2024/8/13 15:38:27.5','2024/07/19 09:36:24']
    # SHIP['SpeedAis'] = [9.58,9.47,8.334,6.94,13,8.82,13.58,13.067,8.1]
    # SHIP['CrossDistDas'] = [3.45,3.1375,3.28,3.07,3.15,2.95,3.075,3.25,3.35]
    # SHIP['Direction'] = ['WS->EN','EN->WS','WS->EN','WS->EN','WS->EN','EN->WS','EN->WS','EN->WS','EN->WS']
    SHIP['MMSI']=SHIP['MMSI'].astype('str')
    SHIP["CrossTimeDas"] = pd.to_datetime(SHIP["CrossTimeDas"],format='%Y/%m/%d %H:%M:%S')
    # SHIP['lat'] = []
    # SHIP['lng'] = []    
    # SHIP['distFromF1'] = [0.4852,0.7085,0.6775,0.8357,0.7691,0.9174,0.8439,0.6879,0.558]

    #AIS原始数据
    PosFile='AIS_DATA_RAW/gulangyu/pos_bj0829_20240715120000_20240816120000_3501.csv'
    StaticFile='AIS_DATA_RAW/gulangyu/static_bj0829_20240715120000_20240816120000_3501.csv'
    ST_UTC8=datetime.datetime.strptime("15/07/24 00:01", "%d/%m/%y %H:%M")
    ET_UTC8=datetime.datetime.strptime("15/08/24 23:59", "%d/%m/%y %H:%M")
        
    dfp=pd.read_csv(PosFile,encoding='GBK')
    dfp['MMSI'] = dfp['MMSI'].astype('str')
    dfp['时间']=pd.to_datetime(dfp['时间'])
    mmsi,Ct = SHIP[['MMSI',"CrossTimeDas"]].iloc[0]
    mmsi_df,lat,lng = getAIS_byCrossTimeDas(mmsi,Ct,dfp)
    
    
    #测试点与点之间的距离是否正确
    for i in range(0,len(SHIP)):
        mmsi,Ct = SHIP[['MMSI',"CrossTimeDas"]].iloc[i]
        _,lat,lng = getAIS_byCrossTimeDas(mmsi,Ct,dfp)
        SHIP['lat'].iloc[i] = lat
        SHIP['lng'].iloc[i] = lng

    SHIP.sort_values(by='CrossDistDas',ascending=True,inplace=True)
    SHIP.reset_index(inplace=True)
    SHIP.dropna(inplace=True)
    # SHIP = SHIP[SHIP['Direction']=='WS->EN']
    # print(SHIP)
    # SHIP = SHIP.iloc[0:-1]
    print(SHIP)
    
    for i in range(1,len(SHIP)):
        distDas = SHIP['CrossDistDas'].iloc[i]-SHIP['CrossDistDas'].iloc[i-1]
        distLat_lng = getDistance(SHIP['lat'].iloc[i-1], SHIP['lng'].iloc[i-1],SHIP['lat'].iloc[i], SHIP['lng'].iloc[i])
        #print(distDas,distLat_lng)
    #fiber0_3d,fiber1_3d = Reposition(SHIP,FIBER_0,FIBER_1)
    fiber0_2d,fiber1_2d,cable_slope,cable_bias = Reposition2D(SHIP,FIBER_0,FIBER_1)
    
    #寻找位于最接近重定位后光缆上方的船舶样本
    print(cable_slope*SHIP.iloc[0]['lat']+cable_bias,SHIP.iloc[0]['lng'])
    SHIP['Reposition_error_lng'] = (cable_slope*SHIP['lat']+cable_bias-SHIP['lng']).abs()
    SHIP.sort_values(by='Reposition_error_lng',ascending=True,inplace=True)
    #print(getDistance(24.450011,118.049537,24.450533,118.049075)*1000)
    fiber0_end_lat_ls=[]
    fiber0_end_lng_ls = []
    fiber1_end_lat_ls=[]
    fiber1_end_lng_ls = []
    
    for i in range(0,len(SHIP)):
        lat0,lng0,dasCh = SHIP.iloc[i]['lat'], SHIP.iloc[i]['lng'],SHIP.iloc[i]['CrossDistDas']
        _, intersect= calculate_distance_and_intersection(cable_slope, cable_bias, lat0,lng0)
        #print('样本在重定位后光缆上的投影经纬度为',intersect)
        Fiber1_end_km = 3.9
        fiber1_end_lat = bisection_method( intersect[0],dasCh,Fiber1_end_km,cable_slope,cable_bias, tol=1e-6)
        fiber1_end_lng = cable_slope*(fiber1_end_lat)+cable_bias
        Fiber0_end_km = 2.73
        fiber0_end_lat = bisection_method( intersect[0],dasCh,Fiber0_end_km,cable_slope,cable_bias, tol=1e-6)
        fiber0_end_lng = cable_slope*(fiber0_end_lat)+cable_bias
        fiber1_end_lat_ls.append(fiber1_end_lat)
        fiber1_end_lng_ls.append(fiber1_end_lng)
        fiber0_end_lat_ls.append(fiber0_end_lat)
        fiber0_end_lng_ls.append(fiber0_end_lng)
    fiber0_end_lat = np.mean(fiber0_end_lat_ls)
    fiber0_end_lng = np.mean(fiber0_end_lng_ls)
    fiber1_end_lat = np.mean(fiber1_end_lat_ls)
    fiber1_end_lng = np.mean(fiber1_end_lng_ls)
    print('重定位后光缆两端经纬度',(fiber0_end_lat,fiber0_end_lng),(fiber1_end_lat,fiber1_end_lng))

    
    #print('重定位后光缆两端经纬度',fiber0_2d,fiber1_2d)
    print('重定位前后光缆两端距离------------------')
    # print(getDistance(fiber0_2d[0], fiber0_2d[1],FIBER_0[0], FIBER_0[1])*1000)
    # print(getDistance(fiber1_2d[0], fiber1_2d[1],FIBER_1[0], FIBER_1[1])*1000)
    
    print('重定位后fiber0端距离',round(getDistance(fiber0_end_lat,fiber0_end_lng,FIBER_0[0], FIBER_0[1])*1000,2))
    print('重定位后fiber1端距离',round(getDistance(fiber1_end_lat,fiber1_end_lng,FIBER_1[0], FIBER_1[1])*1000,2))
    #在map中映射出经纬度

    m = folium.Map(location=[24.4,118.04],width=1200,height=600,zoom_start=12)
    folium.PolyLine([[FIBER_0[0],FIBER_0[1]],[FIBER_1[0],FIBER_1[1]]],color='green').add_to(m)
    folium.PolyLine([[FIBER_0[0],FIBER_0[1]],[FIBER_2[0],FIBER_2[1]]],color='green').add_to(m)
    folium.PolyLine([[FIBER_2[0],FIBER_2[1]],[FIBER_3[0],FIBER_3[1]]],color='green').add_to(m)
    folium.PolyLine([[FIBER_4[0],FIBER_4[1]],[FIBER_3[0],FIBER_3[1]]],color='green').add_to(m)
    folium.PolyLine([[FIBER_5[0],FIBER_5[1]],[FIBER_1[0],FIBER_1[1]]],color='green').add_to(m)

    #读取施工图
    drawing = pd.read_csv('/home/huangwj/DAS/CableReconstruct/Drawing.csv')
    for i in range(0,len(drawing)-1):
        folium.PolyLine([[drawing['Latitude'].iloc[i],drawing['Longtitude'].iloc[i]],[drawing['Latitude'].iloc[i+1],drawing['Longtitude'].iloc[i+1]]],color='yellow').add_to(m)
    #folium.PolyLine([[fiber0_2d[0],fiber0_2d[1]],[fiber1_2d[0],fiber1_2d[1]]],color='orange').add_to(m)
    len_dd = len(drawing)-1
    sd = int(2/3*len_dd)
    ed = int((2/3+0.15)*len_dd)
    for i in range(sd,ed):
        folium.PolyLine([[drawing['Latitude'].iloc[i],drawing['Longtitude'].iloc[i]],[drawing['Latitude'].iloc[i+1],drawing['Longtitude'].iloc[i+1]]],color='black').add_to(m)
    print(drawing['Latitude'].iloc[sd],drawing['Longtitude'].iloc[sd],Fiber0_end_km+getDistance(drawing['Latitude'].iloc[sd],drawing['Longtitude'].iloc[sd],FIBER_0[0], FIBER_0[1]))
    print(drawing['Latitude'].iloc[ed],drawing['Longtitude'].iloc[ed],Fiber0_end_km+getDistance(drawing['Latitude'].iloc[ed],drawing['Longtitude'].iloc[ed],FIBER_0[0], FIBER_0[1]))
    
    folium.PolyLine([[fiber0_end_lat,fiber0_end_lng],[fiber1_end_lat,fiber1_end_lng]],color='orange').add_to(m)

    lat_ls = SHIP['lat'].to_list()
    lng_ls = SHIP['lng'].to_list()
    direct_ls = SHIP['Direction'].to_list()

    for i in range(0,len(lng_ls)):
        if direct_ls[i]=='EN->WS':
            folium.CircleMarker(location = [lat_ls[i], lng_ls[i]],radius=2,color = 'red').add_to(m)
        else:
            folium.CircleMarker(location = [lat_ls[i], lng_ls[i]],radius=2,color = 'blue').add_to(m)

    m.save('map.html')

    #推断光缆走向
    angle = GetCrossAngle(FIBER_3,FIBER_4,FIBER_3,FIBER_2)
    print('施工图中两段光缆的夹角(FIBER_3-FIBER_4-----FIBER_3-FIBER_2):',math.degrees(angle))
    angle = GetCrossAngle(FIBER_2,FIBER_3,FIBER_2,FIBER_0)
    print('施工图中两段光缆的夹角(FIBER_2-FIBER_3-----FIBER_2-FIBER_0):',math.degrees(angle))
    angle = GetCrossAngle(FIBER_0,FIBER_1,FIBER_0,FIBER_2)
    print('施工图中两段光缆的夹角(FIBER_0-FIBER_1-----FIBER_0-FIBER_2):',math.degrees(angle))
    angle = GetCrossAngle(FIBER_1,FIBER_0,FIBER_1,FIBER_5)
    print('施工图中两段光缆的夹角(FIBER_0-FIBER_1-----FIBER_1-FIBER_5):',math.degrees(angle))
    #FIBER_4---FIBER_3---FIBER_2---FIBER_0---FIBER_1---FIBER_5
    print('Fiber4-3 长度',getDistance(FIBER_4[0], FIBER_4[1],FIBER_3[0], FIBER_3[1])*1000)
    print('Fiber3-2 长度',getDistance(FIBER_3[0], FIBER_3[1],FIBER_2[0], FIBER_2[1])*1000)

    print('Fiber2-0 长度',getDistance(FIBER_2[0], FIBER_2[1],FIBER_0[0], FIBER_0[1])*1000)
    print('Fiber0-1 长度',getDistance(FIBER_1[0], FIBER_1[1],FIBER_0[0], FIBER_0[1])*1000)
    print('Fiber1-5 长度',getDistance(FIBER_1[0], FIBER_1[1],FIBER_5[0], FIBER_5[1])*1000)
