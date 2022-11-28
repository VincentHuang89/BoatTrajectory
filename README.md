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
