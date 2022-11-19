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
