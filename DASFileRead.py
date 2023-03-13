
#%%
import os
from datetime import datetime,timedelta

def DasFileRead(StartTime:datetime,EndTime:datetime,Filepath:str):
    T=StartTime
    Times=[]
    result=[]
    #遍历所有从StartTime到EndTime所对应的文件
    while T<=EndTime:  
        filename=T.strftime("%Y%m%d_%H%M")
        files = os.listdir(Filepath)
        for s in files:
            if filename in s:
                Times.append(T)
                result.append(s)
        T=T+timedelta(minutes=1)
    return result,Times



# %%
