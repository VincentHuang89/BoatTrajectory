# BoatTrajectory

20221101

BoatSpeedAngel.py: 用以计算航速和航向范围，输入来自DAS数据的标注斜率

Trajectory_filter.py: 用以描绘DAS数据时空图

WavePatternEval.py: 描绘船行波等相线模式图

20221108

Trajectory_filter：1）对Showdata进行Z-score归一化，再实现imshow，2）调整imshow中的X轴的时间显示格式，删去年月日，保留时分秒。

20221115

Trajectory_filter：1）将关于radon变换以及speed求解的功能重构成RadonWake.py，2）利用Z-score的数学表达式重构了Z-score的功能，3）增加了threshold的功能，在一定范围内的数据置零，以消除海杂波对后续的radon 变换的影响，4）在PlotRadon函数中对radon域的数据进行Z-score归一化，以及实现threshold功能。

20221117

Trajectory_filter: 1)重构PlotDAS功能以及计算DAS图像的SNR, 调整X轴和Y轴的坐标；2)遍历采样率，筛选合适SNR下的speed，输出其平均值。3）修正时间标注功能，将ET调整为ET+minute(1)

20221126

1）对ShowDataSlice进行插值，使空间维度与时间维度一致，提高radon变换的分辨率，2）重构space-time diagram的标注功能，适应于不同的MINCHANNEL和MAXCHANNEL。

20221127

1）ShowDataSlice的空间维度决定了radon域上的角度分辨率。当空间维度较少时（波线不明显，较短时）应在空间维度上通过堆叠矩阵来增加空间维度的信息，因为只需要知道波线的斜率即可，矩阵的堆叠并不会对波线的斜率产生影响，增加了ShowDataSlice的堆叠功能

2）在RadonWake里实现validation函数，以验证所求得波线斜率在ShowDataSlice上与各个波线之间的关系，是否平行。

20221128

拟增加ShowDataSlice的分辨率自动增强功能，因为多个采样率多次采样，其实对角度的提高并没有很好，主要在于欠采样的时候会严重影响波线的斜率，其次，过高的采样率会导致漫长的执行时间。

Trajectory_filter：1）删除过个采样率多次采样的功能，调整为单次采样；2）补充了空间和时间分辨率的自适应增强功能，先通过在时间维度上重采样ShowDataSlice，实现时间维度的自动调整（2000以下），再实现空间维度上的矩阵堆叠，增加空间维度上的信息。

AISData：调整船轨迹过光纤的速度和时间判定，从原来的默认匀速到现在的匀加速。

20230112

Trajectory_filter：按照MMSI来控制输入的数据文件的范围；根据时间控制参数Tstart和Tend来调整PlotDAS函数的时间标注。

20230131:

Paperfigure.py : 用以描绘文章所需的图，包括1）radon 变换前后图片的合并，2）validation of the speed of divergent wave and envelope curve图片的合并，3）simulated divergent wave in DAS，用以比较仿真与真实测量结果

20230311:

WaveVAngel.py: 修改了ROTATE与move函数中对angle，令angle=angle-90，因为wavepattern函数所产生的船行波的方向与光纤垂直。

Trajectory_filter :修改了PlotSimulInDAS函数中的时间标注参数，ST变为ST1等.

20230324:

Paperfigure.py: 重新调整了PlotSimInDAS.png的画图逻辑，Tbias和Wbias依据船行波出现的时刻和位置来自动生成，不需要手动调整，同时调整了Simulated Wake的生成逻辑，仅依靠Wbias来调整尺寸。

Trajectory_filter：增加船行波出现时刻和位置的计算逻辑。

20230611:

AISData.py: 调整光纤部署的起止点，利用百度地图的坐标拾取工具：http://api.map.baidu.com/lbsapi/getpoint/index.html

海图工具：https://www.bilibili.com/video/BV1cT4y1M7kx/?spm_id_from=333.337.search-card.all.click&vd_source=04b165bfde4d781288f5d8beaa41b6ed

![1686489781325](image/README/1686489781325.png)

![1686482650537](image/README/1686482650537.png)

230713: 在AISData.py中，修改getcrossangle函数，将经纬度坐标转成直角坐标系，再利用向量夹角的方式求解船舶路径与光纤的夹角。

230715：在ReVise_Fiber_pos.py中，包含从FTP下载数据的功能、船行波图片的生成以及根据AIS数据进行标注。

230715：在Select_ship.py中，用于根据具体船速，方向等约束条件筛选船只记录，接着进入ReVise_Fiber_pos.py生成相应的图片以及时间标注，用于辅助光纤重定位。

230720：在Trajectory_filter等代码中，修改船只尖峰的定位功能以及数据切片的功能，不再通过t_s等参数，而是预置的数据区域顶点偏置等方式来实现数据切片，从而保存每条船只的数据切片参数。
