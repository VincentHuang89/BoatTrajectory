from tdms_reader_1 import *
import numpy as np
from math import ceil,floor,asin,tan,sin, cos,radians
import pandas as pd
import datetime
import os
from scipy import stats
import matplotlib as mpl
mpl.use('Agg') #必须要写在这两个import中间
import matplotlib.pyplot as plt
from scipy import signal
import math
from WaveVAngel import PlotSimuWaveInDas,FroudeNum,WavePattern,ROTATE,move,cross,WavePattern1

def bandpass_f(some_data,sample_rate,cutoff_f_low,cutoff_f_high,cutoff_order):
    #带通滤波
    data=[]
    cutoff_f0=cutoff_f_low/sample_rate*2 #截止频率10hz
    cutoff_f1=cutoff_f_high/sample_rate*2 
    #cutoff_order=cutoff_order   #配置滤波器的阶数
    if cutoff_f0==0 and cutoff_f1!=0:
        b,a = signal.butter(cutoff_order, cutoff_f1, 'low')
    elif cutoff_f0!=0 and cutoff_f1==0:
        b,a = signal.butter(cutoff_order, cutoff_f0, 'high')
    else:
        b,a = signal.butter(cutoff_order, [cutoff_f0,cutoff_f1], 'bandpass') 
    #分别对所有的通道滤波，滤波后的数据保存在some_data里
    for channel in range(0,some_data.shape[1]):
        data.append(signal.filtfilt(b,a, some_data[:,channel]))
        # if (some_data.shape[1]-channel)%100==0:
        #     print('The rest channels for filtering:{}'.format(some_data.shape[1]-channel))
    return np.array(data).transpose()

def get_sample_png(npy,save_fig_path_name,st,et,dist0_km,dist1_km,Date,fiber1_line,fiber0_line,ch0,tBias,tRange,label = True,line_label = True):
    npy = stats.zscore(npy, axis=0)  #channel-wise z_score
    plt.figure(dpi=1000, figsize=(20, 8))  #dpi=300
    
    # plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = True
    # plt.rcParams['ytick.left'] = plt.rcParams['ytick.labelleft'] = False
    plt.imshow(np.transpose(npy), cmap="bwr", aspect='auto',origin='lower',vmin=-3,vmax=3)
    #plt.imshow(np.transpose(npy), cmap="bwr", aspect='auto',origin='lower',vmin=-np.std(npy)*3,vmax=np.std(npy)*3)
    if label == True:
        plt.grid()
        TimeTicks=3
        xlabel=np.linspace(0,npy.shape[0],TimeTicks+1)
        plt.xticks(xlabel,pd.date_range(st,et+datetime.timedelta(minutes=1),periods=TimeTicks+1).strftime('%H:%M:%S'),rotation = 0,size=20)
        distTicks = 6
        ylabel=np.linspace(0,npy.shape[1],distTicks+1)
        plt.xlabel("Time",fontsize=20)
        plt.yticks(ylabel,np.round(np.linspace(dist0_km,dist1_km,distTicks+1),3),size=20)
        plt.ylabel("Distance(km)",fontsize=20)
        # plt.yticks([])
        #plt.title(Date[:6],size=20)
        if line_label ==True:
            if (fiber1_line<dist1_km)&(fiber1_line>dist0_km):
                fiber1_line_y = (fiber1_line-dist0_km)/(dist1_km-dist0_km)*npy.shape[1]
                plt.axhline(y=fiber1_line_y,color='green',linestyle='--',linewidth=0.5)
            if (fiber0_line<dist1_km)&(fiber0_line>dist0_km):
                fiber0_line_y = (fiber0_line-dist0_km)/(dist1_km-dist0_km)*npy.shape[1]
                plt.axhline(y=fiber0_line_y,color='blue',linestyle='--',linewidth=0.5)
            if (ch0<dist1_km)&(ch0>dist0_km):
                ch0_y = (ch0-dist0_km)/(dist1_km-dist0_km)*npy.shape[1]
                plt.hlines(ch0_y,tBias,tBias+tRange,color='black',linestyle='--',linewidth=1.5)
            
    else:
        plt.axis("off")
    plt.savefig(save_fig_path_name, format='png', bbox_inches='tight', pad_inches=0, dpi=300)  #dpi=300
    plt.close()


def get_sample_and_trace_png(npy,save_fig_path_name,st,et,dist0_km,dist1_km,Date,fiber1_line,fiber0_line,ch0,tBias,tRange,trace,label = True,line_label = True):
    trace_ch = [trace[i][0] for i in range(0,len(trace))]
    trace_t = [trace[i][1] for i in range(0,len(trace))]
    npy = stats.zscore(npy, axis=0)  #channel-wise z_score
    plt.figure(dpi=600, figsize=(10, 6))  #dpi=300
    plt.imshow(np.transpose(npy), cmap="bwr", aspect='auto',origin='lower',vmin=-3,vmax=3)
    if label == True:
        plt.grid()
        TimeTicks=6
        xlabel=np.linspace(0,npy.shape[0],TimeTicks+1)
        plt.xticks(xlabel,pd.date_range(st,et+datetime.timedelta(minutes=1),periods=TimeTicks+1).strftime('%H:%M:%S'),rotation = 30,size=20)
        distTicks = 10
        ylabel=np.linspace(0,npy.shape[1],distTicks+1)
        plt.xlabel("Time",fontsize=20)
        plt.yticks(ylabel,np.round(np.linspace(dist0_km,dist1_km,distTicks+1),3),size=20)
        plt.ylabel("Distance(km)",fontsize=20)
        plt.title(Date[:6],size=15)
        if line_label ==True:
            if (fiber1_line<dist1_km)&(fiber1_line>dist0_km):
                fiber1_line_y = (fiber1_line-dist0_km)/(dist1_km-dist0_km)*npy.shape[1]
                plt.axhline(y=fiber1_line_y,color='green',linestyle='--',linewidth=0.5)
            if (fiber0_line<dist1_km)&(fiber0_line>dist0_km):
                fiber0_line_y = (fiber0_line-dist0_km)/(dist1_km-dist0_km)*npy.shape[1]
                plt.axhline(y=fiber0_line_y,color='blue',linestyle='--',linewidth=0.5)
            if (ch0<dist1_km)&(ch0>dist0_km):
                ch0_y = (ch0-dist0_km)/(dist1_km-dist0_km)*npy.shape[1]
                plt.hlines(ch0_y,tBias,tBias+tRange,color='black',linestyle='--',linewidth=1.5)
            
    else:
        plt.axis("off")
    plt.plot(trace_t,trace_ch,linestyle='-',linewidth='1.5',color = 'black')
    plt.savefig(save_fig_path_name, format='png', bbox_inches='tight', pad_inches=0, dpi=300)  #dpi=300
    plt.close()


def get_sample_png1(npy,save_fig_path_name,st,et,dist0_km,dist1_km,label = True):
    npy = stats.zscore(npy, axis=0)  #channel-wise z_score
    plt.figure(dpi=300, figsize=(20,4))  #dpi=300
    plt.imshow(npy, cmap="bwr", aspect='auto',origin='lower',vmin=-3,vmax=3)
    if label == True:
        plt.grid()
        TimeTicks=6
        ylabel=np.linspace(0,npy.shape[0],TimeTicks+1)
        plt.yticks(ylabel,pd.date_range(st,et+datetime.timedelta(minutes=1),periods=TimeTicks+1).strftime('%H:%M:%S'),rotation = 30,size=10)
        distTicks = 17
        plt.ylabel("Time (2024/08/14)",fontsize=15)
        xlabel=np.linspace(0,npy.shape[1],distTicks+1)
        plt.xticks(xlabel,np.round(np.linspace(dist0_km,dist1_km,distTicks+1),3),size=15)
        plt.xlabel("Distance(km)",fontsize=15)
    else:
        plt.axis("off")
    plt.savefig(save_fig_path_name, format='png', bbox_inches='tight', pad_inches=0, dpi=300)  #dpi=300
    plt.close()
    

def get_data_from_submarineFbier(npy,dist0_km,dist1_km,Channel_spacing,zero_offset):
    #npy.shape[1]为channel
    ch0= floor((dist0_km*1000-zero_offset)/Channel_spacing)
    ch1= ceil((dist1_km*1000-zero_offset)/Channel_spacing)
    return npy[:,ch0:ch1]
    
def get_time_from_fileName(filename):
    return datetime.datetime.strptime(filename[-24:-9],'%Y%m%d_%H%M%S')+datetime.timedelta(hours=8)



def get_one_ch(data,ch,tBias,tRange,dist0_km,dist1_km,incr_ch=0,SamplesPerSec=100,time_shift=0):
    data = stats.zscore(data, axis=0)
    ch_idx = int((ch-dist0_km)/(dist1_km-dist0_km)*data.shape[1])-incr_ch
    st_idx = int(tBias*SamplesPerSec)-time_shift
    et_idx = st_idx+int((tRange)*SamplesPerSec)
    return data[st_idx:et_idx,ch_idx],(ch_idx,st_idx,et_idx)



def get_one_ch_idx(data,ch,st_idx,et_idx,dist0_km,dist1_km,incr_ch=0):
    data = stats.zscore(data, axis=0)
    ch_idx = int((ch-dist0_km)/(dist1_km-dist0_km)*data.shape[1])-incr_ch
    return data[st_idx:et_idx,ch_idx],(ch_idx,st_idx,et_idx)

    
def get_chs_png(data_chs,save_fig_path_name):
    plt.figure(dpi=300, figsize=(20,4))  #dpi=300
    for i in range(0,len(data_chs)):
        if i ==0:
            plt.plot(data_chs[i]-i*2,linestyle = '--')
        else:
            plt.plot(data_chs[i]-i*2)
    plt.savefig(save_fig_path_name, format='png', bbox_inches='tight', pad_inches=0, dpi=300)
 
def get_some_chs(data,ch,tBias,tRange,dist0_km,dist1_km,decr_chs=10,SamplesPerSec=100 ):
    datas = []
    for decr_ch in np.arange(0,decr_chs):
        datas.append(get_one_ch(data,ch,tBias,tRange,dist0_km,dist1_km,decr_ch,SamplesPerSec)[0])
    get_chs_png(datas,save_fig_path_name)
    return datas

def track_some_chs(data,ch,tBias,tRange,dist0_km,dist1_km,time_shifts,decr_chs,SamplesPerSec,save_fig_path_name ):
    datas = []
    time_shifts = [0]+time_shifts
    time_shifts = [sum(time_shifts[:i+1]) for i in range(len(time_shifts))]
    for decr_ch,time_shift in zip(np.arange(0,decr_chs),time_shifts):
        datas.append(get_one_ch(data,ch,tBias,tRange,dist0_km,dist1_km,decr_ch,SamplesPerSec,time_shift)[0])
    get_chs_png(datas,save_fig_path_name)
    return datas




def track_some_chwise(data,ch,tBias,tRange,dist0_km,dist1_km,decr_chs,SamplesPerSec,save_fig_path_name ):
    datas = []
    delta_sec = 1
    st = []
    trace = []
    for decr_ch in np.arange(0,decr_chs):
        if decr_ch==0:
            time_shift = 0
            data_x,mes = get_one_ch(data,ch,tBias,tRange,dist0_km,dist1_km,decr_ch,SamplesPerSec,time_shift)      
            data_len = len(data_x)      
            datas.append(data_x)
            st.append(mes[1])
            trace.append((mes[0],(mes[1]+mes[2])//2))
        else:
            data_y,mes = get_one_ch_idx(data,ch,mes[1]-delta_sec*SamplesPerSec,mes[2]+delta_sec*SamplesPerSec,dist0_km,dist1_km,decr_ch)
            #分析互相关
            time_shift = get_delay(data_x,data_y)
            data_x,mes = get_one_ch_idx(data,ch,mes[1]+time_shift,mes[1]+time_shift+data_len,dist0_km,dist1_km,decr_ch)
            datas.append(data_x)
            st.append(mes[1])
            trace.append((mes[0],(mes[1]+mes[2])//2))

    get_chs_png(datas,save_fig_path_name)
    time_shifts = [st[i+1]-st[i] for i in range(0,len(st)-1)]
    return datas,time_shifts,trace        
        
    

 
def get_some_chs_shift(data_chs,time_shifts,save_fig_path_name):
    
    time_shifts = [0]+time_shifts
    
    time_shifts = [sum(time_shifts[:i+1]) for i in range(len(time_shifts))]
    plt.figure(dpi=300, figsize=(20,4))  #dpi=300
    for i in range(0,len(data_chs)):
        if i ==0:
            plt.plot(np.arange(time_shifts[i],time_shifts[i]+len(data_chs[i]),1),data_chs[i]-i*2,linestyle = '--')
        else:
            plt.plot(np.arange(time_shifts[i],time_shifts[i]+len(data_chs[i]),1),data_chs[i]-i*2)
    plt.savefig(save_fig_path_name, format='png', bbox_inches='tight', pad_inches=0, dpi=300)


def get_delay(data_x,data_y):
    #data_x: 短序列，data_y:长序列
    if len(data_x)==len(data_y):
        correlation = np.correlate(data_x,data_y, mode='full')
        max_corr_index = np.argmax(correlation)
        time_shift = max_corr_index - len(data_x) + 1
    else:
        correlation = np.correlate(data_y,data_x, mode='valid')
        max_corr_index = np.argmax(correlation)
        time_shift = max_corr_index
    return time_shift

def get_delays(data_chs):
    time_shifts= []
    for i in range(0,len(data_chs)-1):
        data_x = data_chs[i]
        data_y = data_chs[i+1]
        time_shifts.append(get_delay(data_x,data_y))
    return time_shifts

def PLotSimWaveCrest(UpperBound,LowerBound,angle,v,T,N,a,delta_T):
    crossp1=[]
    crossp2=[]
    t_start=0
    t_start1=0
    for t in np.arange(0,T,delta_T):
        dist=v*t #移动的距离 dist=v*t,除以矫正系数
        X,Y,Y1,ALPHA=WavePattern1(a,N,0.4,2,0.01)   #0.5
        #Y1=list(np.array(Y1)-0)   #调整单边散波波形的位置
        alpha=max(ALPHA)
        Attenuation0 = sin(radians(angle - (alpha)))**2   #乘以一定系数
        Attenuation1 = sin(radians(180-angle - (alpha)))**2
        RX,RY=ROTATE(angle,X,Y)
        RX1,RY1=ROTATE(angle,X,Y1)
        RMX,RMY=move(RX,RY,angle,dist)
        RMX1,RMY1=move(RX1,RY1,angle,dist)
        cp=cross(RMX,RMY)
        #Apply the upper and the lower bound of the wave crest to fit the DAS data
        if len(cp)>0:
            if cp[0]>UpperBound or cp[0]<LowerBound: 
                cp=[]
        if len(cp) and t_start==0:
            t_start=t
        crossp1=crossp1+cp

        cp=cross(RMX1,RMY1)
        if len(cp)>0:
            if cp[0]>UpperBound or cp[0]<LowerBound: 
                cp=[]
        if cp and t_start1==0:
            t_start1=t
        crossp2=crossp2+cp
        x1=np.arange(0,len(crossp1),1)
        x2=np.arange(0,len(crossp2),1)
        x1=x1*delta_T+t_start
        x2=x2*delta_T+t_start1
    return x1,x2,crossp1,crossp2,Attenuation0,Attenuation1


def PlotSimulInDAS(REGION,DownSampleRate,v,h,angle,A,ShowData,ST,ET,MINCHANNEL,channel_spacing,WLen_Scale,Wbias,Tbias,UpperBound,LowerBound,SHIP):
    '''
    Depict the simulated ship wake in the DAS figure to analyze the difference between the simulated  and measured divergent wave
    '''
    Tim,Pos = REGION['Wave_peak']
    T=60
    N=FroudeNum(v,h)
    print('Froude Number is ',N)

    plt.figure(dpi=400,figsize=(15,8))
    
    delta_T=0.1
    for a in np.arange(1,A+1,4):
        x1,x2,crossp1,crossp2,Attenuation0,Attenuation1=PLotSimWaveCrest(UpperBound,LowerBound,angle,v,T,N,a,delta_T)
        #print("散波速度",(crossp1[-1]-crossp1[-2])/(x1[-1]-x1[-2]))
        T_bias = (a-1)*1  #控制各条波线之间的间距
        #Scale the wave length because the WavePattern function returns the dimensionless wake pattern
        x1=WLen_Scale*x1
        x2=WLen_Scale*x2
        crossp1=WLen_Scale*(np.array(crossp1))
        crossp2=WLen_Scale*(np.array(crossp2))

        #Channel sampling and sampling of the simulated wave crest
        crossp1=crossp1/channel_spacing+Pos
        crossp2=crossp2/channel_spacing+Pos
        x1=x1*DownSampleRate+Tim
        x2=x2*DownSampleRate+Tim
       
        plt.xlim((0,ShowData.shape[0]))
        plt.ylim(((0,ShowData.shape[1])))
        plt.plot(x1,crossp1,lw=1.5,linestyle='--',alpha=Attenuation0,color='lime')
        plt.plot(x2,crossp2,lw=1.5,linestyle='--',alpha=Attenuation1,color='lime')

    
    ShipWave_top_pos = REGION['Wave_peak']
    #ax1 = fig1.add_subplot(1,1,1) 
    plt.imshow(np.transpose(ShowData), cmap="bwr", aspect='auto',origin='lower',vmin=-3,vmax=3,extent=[0, ShowData.shape[0], 0, ShowData.shape[1]]) # cmap=''bwr,
    plt.colorbar()

    #plt.scatter(ShipWave_top_pos[0], ShipWave_top_pos[1], color='yellow', marker='o')
    #计算坐标轴的刻度大小,合计10个时间戳(计算有问题，需要考虑数据的实际距离以及截断)
    #根据MINCHANNEL和MAXCHANNEL重置Y轴
    TimeTicks = 7
    xlabel = np.linspace(0, ShowData.shape[0] - 1, TimeTicks)
    xtick_labels = pd.date_range(ST.strftime("%Y%m%d %H%M%S"), ET.strftime("%Y%m%d %H%M%S"), periods=TimeTicks).strftime('%H:%M:%S')
    TimeTicks=7
    xlabel=np.linspace(0,ShowData.shape[0]-1,TimeTicks)
    plt.xticks(xlabel,xtick_labels,rotation = 0,size=15)
    ylabel=np.linspace(0,ShowData.shape[1],5)
    plt.xlabel("Time",fontsize=15)
    plt.yticks(ylabel,np.round((ylabel)*channel_spacing/1000+MINCHANNEL,2),size=15)
    plt.ylabel("Distance(Km)",fontsize=15)

    plt.savefig('PlotSimulInDAS.png',bbox_inches='tight')
    plt.savefig(f'/home/huangwj/DAS/BoatTrajectory/Paperfig/PlotSimulInDAS_{SHIP}.pdf',bbox_inches='tight')



        
        
    
if __name__=='__main__':

    
    Fi = 6
    DataPath = 'DASforAIS_gly' #w文件存放路径
    SubFolder = ['2408111042',"2407221131",'2408080729','2407190150','2408021250','2408041728','2408140046','2408061207',  '2407201154','2407271808','2408101544','2408091028','2407301211','2407181529','2408050933','2408061759','2408131538','2407191920','2407241012','2407301907','2408150837','2407190936','2408150832','2408140002','2408051246','2407171009','2408071008','2407190950','2408010804','2408130340','2407301044'] #
    TdmsPath = os.path.join(DataPath,SubFolder[Fi])
    Fileset = os.listdir(TdmsPath)
    Files_df = pd.DataFrame(Fileset,columns=['filename'])
    Files_df.sort_values(by='filename',ascending=True,inplace=True)
    Files_df['Time'] = Files_df.apply(lambda x:get_time_from_fileName(x['filename']),axis=1)
    print(Files_df)
    FileSet = Files_df.iloc[13:17] #len(Files_df)
    st = FileSet['Time'].iloc[0]
    et = FileSet['Time'].iloc[-1]
    print(FileSet,st,et)
    #目标海缆的起止位置，相对A端机房
    dist0_km=0
    dist1_km=14.825 #14.825
    dist0_km=2.3
    dist1_km=4
    
    ch0_km =3.02
    tBias = 145 #57seconds
    tRange = 10
    
    Fiber1_line = 3.026  #按照施工图估算出来
    Fiber0_line = 2.848
    
    
    all_data_downsampled = pd.DataFrame()
    SamplesPerSec = 100
    for i in range(0,len(FileSet)):
        FilePath= os.path.join(TdmsPath,FileSet['filename'].iloc[i])

        if i==0:
            tdms = TdmsReader(FilePath)
            props = tdms.get_properties()
            zero_offset = props.get('Zero Offset (m)')
            channel_spacing = props.get('SpatialResolution[m]') * props.get('Fibre Length Multiplier')
            n_channels = tdms.fileinfo['n_channels']
            depth = zero_offset + np.arange(n_channels) * channel_spacing
            fs = props.get('SamplingFrequency[Hz]')
            print(f'Zero_offset is:{zero_offset}')
            print('Channel_spacing is:{:.4f}m'.format(channel_spacing))
            print('Number of channels in file: {0}'.format(n_channels))
            print('Time samples in file: {0}'.format(tdms.channel_length))
            print('Sampling frequency (Hz): {0}'.format(fs))
            print(f'Total length: {np.max(depth)}m')
            down_sample_factor = int(props.get('SamplingFrequency[Hz]')/SamplesPerSec)
            first_channel = 0  
            last_channel = n_channels - 1  
            first_time_sample = 00 
            last_time_sample = tdms.channel_length - 1 
            #some_data = tdms.get_data(first_channel, last_channel, first_time_sample, last_time_sample)   
            some_data = pd.DataFrame(tdms.get_data(first_channel, last_channel, first_time_sample, last_time_sample))
            #some_data = pd.DataFrame(bandpass_f(some_data.to_numpy(),SamplesPerSec,0,1,3))
            data_downsampled = some_data.groupby(some_data.index // down_sample_factor).mean()  # 平均采样
            del some_data  # 释放内存
            all_data_downsampled = pd.concat([all_data_downsampled, data_downsampled], axis=0)
        else:
            tdms = TdmsReader(FilePath)
            some_data = pd.DataFrame(tdms.get_data(first_channel, last_channel, first_time_sample, last_time_sample)) 
            #some_data = pd.DataFrame(bandpass_f(some_data.to_numpy(),SamplesPerSec,0,1,3))
    
       
            data_downsampled = some_data.groupby(some_data.index // down_sample_factor).mean()  # 平均采样
            del some_data  # 释放内存
            all_data_downsampled = pd.concat([all_data_downsampled, data_downsampled], axis=0)
        print(FilePath)
    all_data_np = all_data_downsampled.to_numpy()
    
    data_np = get_data_from_submarineFbier(all_data_np,dist0_km,dist1_km,channel_spacing,zero_offset)
    print(data_np.shape)
    save_fig_path_name = f'gly_fig/{SubFolder[Fi]}.png'
    #get_sample_png(data_np,save_fig_path_name,st,et,dist0_km,dist1_km,SubFolder[Fi],Fiber1_line,Fiber0_line,ch0_km,tBias*SamplesPerSec,tRange*SamplesPerSec,True,False)
    get_sample_png1(data_np,save_fig_path_name,st,et,dist0_km,dist1_km,True)

    save_fig_path_name = f'gly_fig/{SubFolder[Fi]}_ch.png'
    decr_chs = 10
    data_chs = get_some_chs(data_np,ch0_km,tBias,tRange,dist0_km,dist1_km,decr_chs,SamplesPerSec=100 )
    #分析连续通道间的时间差异
    #time_shift = get_delay(data_chs[0],data_chs[1])
    time_shifts = get_delays(data_chs)  #依照不同通道选取数据
    
    save_fig_path_name = f'gly_fig/{SubFolder[Fi]}_ch_shift.png'
    get_some_chs_shift(data_chs,time_shifts,save_fig_path_name)
    
    #按照上述time_shifts来重新获取各个通道的数据，跟踪波峰
    save_fig_path_name = f'gly_fig/{SubFolder[Fi]}_ch_track.png'

    data_chs_shift= track_some_chs(data_np,ch0_km,tBias,tRange,dist0_km,dist1_km,time_shifts,decr_chs,SamplesPerSec,save_fig_path_name )
    time_shifts = get_delays(data_chs_shift)  #依照不同通道选取数据
    decr_chs=40
    save_fig_path_name = f'gly_fig/{SubFolder[Fi]}_ch_shifts.png'
    datas,time_shifts,trace = track_some_chwise(data_np,ch0_km,tBias,tRange,dist0_km,dist1_km,decr_chs,SamplesPerSec,save_fig_path_name )
    per_time_shifts = [round(time_shifts[i+1]/time_shifts[i]-1,1) for i in range(0,len(time_shifts)-1)]
    save_fig_path_name = f'gly_fig/{SubFolder[Fi]}.png'
    #get_sample_and_trace_png(data_np,save_fig_path_name,st,et,dist0_km,dist1_km,SubFolder[Fi],Fiber1_line,Fiber0_line,ch0_km,tBias*SamplesPerSec,tRange*SamplesPerSec,trace,True,True)
    speed_proj = [round(SamplesPerSec*channel_spacing/(trace[i+1][1]-trace[i][1]),1) for i in range(0,len(trace)-1)]
    print(speed_proj)
    #2407180950 的船速为9.58
    V_ship = 8.334
    fai = [asin(min(1,V_ship/speed_proj[i])) for i in range(0,len(speed_proj))]
    print([math.degrees(f) for f in fai])
    angle_change= [math.degrees(f)-90 for f in fai]
    print(angle_change)
    change_unit_x = [channel_spacing*cos(math.radians(a)) for a in angle_change]
    change_unit_y = [channel_spacing*sin(math.radians(a)) for a in angle_change]
    route_x = [sum(change_unit_x[0:i]) for i in range(0,len(change_unit_x))]
    route_y = [sum(change_unit_y[0:i]) for i in range(0,len(change_unit_y))]
    save_fig_path_name = f'gly_fig/{SubFolder[Fi]}_route.png'
    plt.figure(dpi=300, figsize=(20,4))  #dpi=300
    plt.plot(route_x,route_y)
    plt.savefig(save_fig_path_name, format='png', bbox_inches='tight', pad_inches=0, dpi=300)
