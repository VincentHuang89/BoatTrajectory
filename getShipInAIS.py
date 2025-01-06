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

ST_UTC8=datetime.datetime.strptime("15/07/24 00:01", "%d/%m/%y %H:%M")
ET_UTC8=datetime.datetime.strptime("15/08/24 23:59", "%d/%m/%y %H:%M")
filename='AIS_SHIP_DATA/FiberBoatMessage_gly_240715_240815.csv'

SHIP= pd.read_csv(filename,index_col=0)
SHIP['CrossTime'] = pd.to_datetime(SHIP['CrossTime'])

SHIP_filter = SHIP[(SHIP['Time_0_1(min)']<=1)&(SHIP['err_speed'].abs()<2)&(SHIP['CrossSpeed']>=6)]  #
#SHIP_filter = SHIP_filter[(SHIP_filter['Angle']<90)&(SHIP_filter['disFromEnd/Km']>0.65)]
SHIP_filter.sort_values(by=['CrossSpeed','disFromEnd/Km'],ascending=False,inplace=True)
print(SHIP_filter.columns)
SHIP_filter[['MMSI', 'CrossTime','CrossSpeed', 'lat', 'lng','disFromEnd/Km', 'tra_direction', 'Angle','ShipName',  'length', 'width', 'depth', 'type',]].to_csv('AIS_SHIP_DATA/FiberBoatMessage_gly_filter.csv')
print(SHIP_filter[['MMSI', 'CrossTime','CrossSpeed', 'lat', 'lng','disFromEnd/Km', 'tra_direction', 'Angle','ShipName',  'length', 'width', 'depth', 'type',]])
