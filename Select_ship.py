import pandas as pd
import os
os.chdir('/home/huangwj/DAS/BoatTrajectory/')

flag=2
if flag ==1:
#2022 7月15到8月15 AIS数据
    filename='AIS_SHIP_DATA/FiberBoatMessage_220715_220815.csv'
elif flag ==2:
#2021 9月份AIS数据
    filename='AIS_SHIP_DATA/FiberBoatMessage_210901_210930.csv'

elif flag ==3:
#2021 10-01到12-01 的AIS数据
    filename='AIS_SHIP_DATA/FiberBoatMessage_211001_211201.csv'

elif flag ==4:
#2022 04-01  07-01 的AIS数据
    filename='AIS_SHIP_DATA/FiberBoatMessage_220401_220701.csv'

else:
    filename='AIS_SHIP_DATA/FiberBoatMessage_220901_230501.csv'

SHIP_DATA = pd.read_csv(filename,index_col=0)
print(SHIP_DATA.columns)


#筛选速度大于11m/s 以及位置小于11，角度在80-100左右的船只

SHIP_DATA_selected = SHIP_DATA[(SHIP_DATA['CrossSpeed']>10)&(SHIP_DATA['disFromEnd/Km']>4)] #&(SHIP_DATA['Angle']>70)&(SHIP_DATA['Angle']<120)
SHIP_DATA_selected.to_csv('AIS_SHIP_DATA/SHIP_selected.csv')
print(SHIP_DATA_selected)