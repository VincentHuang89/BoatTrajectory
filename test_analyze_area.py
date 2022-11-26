import numpy as np
import matplotlib.pyplot as plt
from math import pi, e
from tdms_reader_1 import *
import scipy.signal as ss
import glob

def FK_spectrum_by_2DFFT(data, dt, dx, r_type = 'full'):
    """计算F-K谱
    :param data: 二维地震数据(t-d)
    :param dt: 时间采样间隔
    :param dx: channel间距
    :param L: 平滑窗口长度(暂时没用)
    :param r_type: 选择输出全谱/单边谱
    :return: S(F-K谱),f(频率范围),k(波数范围)
    """
    data = np.array(data)
    [nt, nx] = data.shape

    # 计算nk,nf，FFT补零
    i = 0
    while (2 ** i) <= nx:
        i = i + 1
    nk = 4 * 2 ** i
    j = 0
    while (2 ** j) <= nt:
        j = j + 1
    nf = 4 * 2 ** j
    nf = nt
    nk = nx

    # 2D FFT
    S = np.fft.fftshift(np.fft.fft2(data, (nf, nk)))

    # 设置汉明窗，减少旁瓣泄露
    # H1 = np.hamming(L)
    # H = H1.reshape(L, -1) * H1.reshape(1, L)
    # S = ss.convolve2d(S, H, boundary='symm', mode='same')  # 汉明平滑

    if r_type == 'full':
        f = np.arange(-nf / 2, nf / 2, 1)
        k = np.arange(-nk / 2, nk / 2, 1)
    elif r_type == 'half':
        S = S[nf // 2:nf, ::-1]  # 保留一二象限并左右翻转
        f = np.arange(0, nf / 2, 1)
        k = np.arange(-nk / 2, nk / 2, 1)
    else:
        print('Invalid return type！')

    # 归一化
    S = S / nf / nk
    f = f / nf / dt
    k = k / nk / dx
    return S, f, k


def Dt_by_2DiFFT(FK_S, f_max, k_max, r_type='full'):
    """通过F-K谱，计算Distance-time图
    :param FK_S: 二维地震数据(F-K谱)
    :param f_max: 最大频率
    :param k_max: 最大波数
    :param r_type: 提示输入全谱/单边谱
    :return: Dt_data: 二维地震数据(Distance-time图)
    :return dt: 时间采样间隔
    :return dx: 空间采样间距
    """
    dt = 1/f_max/2
    dx = 1/k_max/2

    FK_S = np.array(FK_S)
    [nf, nk] = FK_S.shape

    if r_type == 'full':
        pass
    elif r_type == 'half':
        #  todo:利用对称性翻转补全F-K谱
        pass
    else:
        print('Invalid return type！')

    # 2D iFFT
    Dt_data = np.fft.ifft2(np.fft.ifftshift(FK_S), (nf, nk))

    return Dt_data, dt, dx


if __name__ == '__main__':
    # print("TDMS Reader demo.")
    # 设定工作路径
    path = r'/home/huangwj/DAS/BoatTrajectory'
    os.chdir(path)
    print('当前工作路径是：' + os.getcwd())  # 显示当前路径
    # 设定提取数据文件
    file_path = 'DataforAIS/sanjiao-1k-4m_UTC_20220724_130942.090.tdms'

    print('File: {0}'.format(file_path))

    tdms = TdmsReader(file_path)

    props = tdms.get_properties()

    zero_offset = props.get('Zero Offset (m)')

    channel_spacing = props.get('SpatialResolution[m]') * props.get('Fibre Length Multiplier')
    n_channels = tdms.fileinfo['n_channels']
    depth = zero_offset + np.arange(n_channels) * channel_spacing
    fs = props.get('SamplingFrequency[Hz]')

    print('Number of channels in file: {0}'.format(n_channels))
    print('Time samples in file: {0}'.format(tdms.channel_length))
    print('Sampling frequency (Hz): {0}'.format(fs))

    """ 选择数据范围 """
    first_channel = 2600  # 选择开始的位置
    last_channel = 2820-1  # 选择结束位置
    first_time_sample = 40*1000  # 选择开始时间，单位ms
    last_time_sample = 45*1000 - 1  # 选择截止时间，单位ms
    filter_range = '50-70'  # Hz

    some_data = tdms.get_data(first_channel, last_channel, first_time_sample, last_time_sample)
    print('Size of data loaded: {0}'.format(some_data.shape))

    plot_derate = 1  # 绘制的降采率
    plot_data = some_data[0::plot_derate, :]

    """画原始数据"""
    plt.clf()
    plt.figure(dpi=100)
    vmin1 = -3
    vmax1 = 3
    from scipy import stats
    plot_data = stats.zscore(plot_data, axis=None)
    # vmin = -np.std(plot_data) * 1
    # vmax = np.std(plot_data) * 5

    extent = [first_time_sample / fs, first_time_sample / fs+(plot_data[::, 0].shape[0] - 1) * plot_derate / fs,
              depth[first_channel], depth[last_channel]]
    img1 = plt.imshow(plot_data.T, cmap='bwr', aspect='auto', vmin=vmin1, vmax=vmax1, extent=extent,
                      origin='lower')
    plt.xlabel('Time (s)', fontsize=12)
    plt.ylabel('Distance (m)', fontsize=12)

    os.chdir(path + r'/F-K谱')
    plt.colorbar(img1)
    plt.title(file_path)

    view_files = file_path
    time_offset = (last_time_sample-first_time_sample+1)//1000  # 单位：s
    date_time = view_files.split("_")
    date_time = date_time[2] + date_time[3][:6]
    output_name = f'view_{date_time}_{int(time_offset)}s_DS{plot_derate}_ch{first_channel}-{last_channel}-{round(channel_spacing, 2)}'
    output_number = len(glob.glob(f'{output_name}*'))
    plt.savefig(f'{output_name}_{output_number}.png', format='png', bbox_inches='tight', pad_inches=0.1)
    plt.clf()

    """ 画F-K谱并保存 """
    plt.figure(dpi=100)

    [FK_S, f, k] = FK_spectrum_by_2DFFT(plot_data, plot_derate / fs, channel_spacing)

    # cmap为图片显示的颜色，autumn:红-橙-黄 gray:黑-白 hot:黑-红-黄-白 jet:蓝-青-黄-红(prefer)
    plot_zoom = False
    if plot_zoom:  # 放大到[-fm,fm]Hz,f区域按比例scale_k中心对称缩放
        [nf, nk] = FK_S.shape
        scale_k = 1/2  # 0.25/5  # 缩放k区间比例
        fm = 5  # Hz
        idx_fm = int(nf / fs * plot_derate*fm)

        if scale_k<1:
            FK_S_plot = FK_S[int(nf/2)-idx_fm:int(nf/2)+idx_fm+1:, int(nk*(1-scale_k)/2):-int(nk*(1-scale_k)/2)+1]
            K_list = k[int(nk * (1-scale_k)/2):-int(nk * (1-scale_k)/2)+1]
        elif scale_k==1:
            FK_S_plot = FK_S[int(nf/2)-idx_fm:int(nf/2)+idx_fm+1:, :]
            K_list = k

        FK_S_plot = stats.zscore(np.abs(FK_S_plot), axis=None)
        vmin = -3
        vmax = 3

        F_list = f[int(nf/2)-idx_fm:idx_fm+int(nf/2)+1]

        extent2 = [F_list[0], F_list[-1], K_list[0], K_list[-1]]  # k-f谱
        print(extent)
        # img1 = plt.imshow(FK_S_plot.T, cmap='bwr', aspect='auto', extent=extent, vmin=vmin, vmax=vmax, origin='lower')

    else:
        extent2 = [f[0], f[-1], k[0], k[-1]]  # k-f谱
        print(extent)
        vmax = 3
        F_list = f
        K_list = k
        FK_S_plot = stats.zscore(np.abs(FK_S), axis=None)

    img1 = plt.imshow(FK_S_plot.T, cmap='bwr', aspect='auto', extent=extent2, vmin=-vmax, vmax=vmax, origin='lower')

    plt.xlabel('Frequency (Hz)', fontsize=12)
    plt.xlim([0, 5])
    plt.ylabel('Wave number (1/m)', fontsize=12)

    plt.colorbar(img1)
    plt.title('F-K Spectrum')

    # 画色散曲线
    line_style = '--'
    H = 16.6  # 水深(m)
    F_mesh, K_mesh = np.meshgrid(F_list, K_list)
    z = 2 * pi * np.power(F_mesh, 2) - 9.8 * K_mesh * np.tanh(2 * pi * K_mesh * H)
    # pc2 = plt.contour(F_mesh, K_mesh, z, 0, colors='w', linestyles=line_style, linewidths=1)
    # plt.plot([], color='w', linestyle=line_style, label=f'H={H}')
    # plt.legend()

    # plt.xlim([-0.5, 0.5])

    # H = 8  # 水深(m)
    # z = 2 * pi * np.power(F_mesh, 2) - 9.8 * K_mesh * np.tanh(2 * pi * K_mesh * H)
    # pc2 = plt.contour(F_mesh, K_mesh, z, 0, colors='r', linestyles=line_style, linewidths=1)
    # plt.plot([], color='r', linestyle=line_style, label=f'H={H}')
    #
    # H = 9  # 水深(m)
    # z = 2 * pi * np.power(F_mesh, 2) - 9.8 * K_mesh * np.tanh(2 * pi * K_mesh * H)
    # pc2 = plt.contour(F_mesh, K_mesh, z, 0, colors='k', linestyles=line_style, linewidths=1)
    # plt.plot([], color='k', linestyle=line_style, label=f'H={H}')


    plt.savefig(f'FK-{output_name}_{output_number}.png', format='png', bbox_inches='tight', pad_inches=0.1)
    plt.clf()

    """画C_p-K图"""
    # plt.figure()
    # line_style = '-'
    # # print(K_list[int(nk * scale_k/2)+
    # sel_1 = True
    # if sel_1:  # 取第一象限
    #     F_mesh, K_mesh = np.meshgrid(F_list[idx_fm:], K_list[int(nk * scale_k/2)+1:])  # 取第一象限
    #     C_p_mesh = F_mesh/K_mesh
    #     levels = np.linspace(0, 10, 4000)
    #     extent = [F_list[idx_fm], F_list[-1], np.min(C_p_mesh), np.max(C_p_mesh)]
    #     # img1 = plt.imshow(FK_S_plot[idx_fm:, int(nk * scale_k/2)+1:].T, aspect='auto', extent=extent, vmin=vmin, vmax=vmax, origin='lower')
    #     pc = plt.contour(F_mesh, C_p_mesh, FK_S_plot[idx_fm:, int(nk * scale_k/2)+1:].T, levels=levels, linestyles=line_style)
    #     plt.colorbar()
    #
    #     plt.xlim([0, 0.5])
    #     plt.ylim([0, 20])
    #
    # else:  # 取第二象限
    #     F_mesh, K_mesh = np.meshgrid(F_list[:idx_fm+1], K_list[int(nk * scale_k / 2) + 1:])  # 取第二象限
    #     C_p_mesh = F_mesh / K_mesh
    #     levels = np.linspace(0, 10, 4000)
    #     extent = [F_list[0], F_list[idx_fm], np.min(C_p_mesh), np.max(C_p_mesh)]
    #     pc = plt.contour(F_mesh, C_p_mesh, FK_S_plot[:idx_fm+1, int(nk * scale_k / 2) + 1:].T, levels=levels,
    #                      linestyles=line_style)
    #     plt.colorbar()
    #
    #     plt.xlim([-0.5, 0])
    #     plt.ylim([-20, 0])
    #
    # z = 2 * pi * F_mesh * C_p_mesh - 9.8 * np.tanh(2 * pi * H * F_mesh / C_p_mesh)
    # # 作图
    # # 绘制等高线,8表示等高线数量加1
    # pc2 = plt.contour(F_mesh, C_p_mesh, z, 0, colors='r', linestyles=line_style)
    # plt.plot([], color='r', linestyle=line_style, label=f'H={H}')
    #
    # plt.xlabel('Frequency (Hz)', fontsize=12)
    # plt.ylabel('Phase velocity (m/s)', fontsize=12)
    #
    # # plt.ylim([0, 2000])
    #
    # plt.title('Cp-F Spectrum')
    # # plt.clabel(pc, inline=True, fontsize=6)  # 添加等高线z标签
    #
    # plt.legend()
    #
    # plt.savefig(f'Cp-F-{output_name}_{output_number}.png', format='png', bbox_inches='tight', pad_inches=0.1)


    '''IFFT'''
    [nf, nk] = FK_S.shape

    if filter_range.split('-')[0]=="0":  # 'lowpass'
        mask = np.zeros((nf, nk), np.uint8)
        scale_k = 1  # 3 / 4  # 0.25/5  # 缩放k区间比例
        if filter_range.split('-')[1]!="m":
            fm = float(filter_range.split('-')[1])  # Hz
        else:
            fm = f[-1]
        idx_fm = float(nf / fs * plot_derate * fm)
        if scale_k < 1:
            mask[int(nf / 2) - idx_fm:int(nf / 2) + idx_fm + 1:,
            int(nk * (1 - scale_k) / 2):-int(nk * (1 - scale_k) / 2) + 1] = 1
        elif scale_k == 1:
            mask[int(nf / 2- idx_fm) :int(nf / 2+ idx_fm)+1:, :] = 1


        fmin = 0
        fmax = fm
    elif filter_range[-1]=="m":  # 'highpass'
        mask = np.ones((nf, nk), np.uint8)
        scale_k = 1  # 3 / 4  # 0.25/5  # 缩放k区间比例
        fm = float(filter_range.split('-')[0])  # Hz
        idx_fm = int(nf / fs * plot_derate * fm)
        if scale_k < 1:
            mask[int(nf / 2) - idx_fm:int(nf / 2) + idx_fm + 1:,
            int(nk * (1 - scale_k) / 2):-int(nk * (1 - scale_k) / 2) + 1] = 0
        elif scale_k == 1:
            mask[int(nf / 2) - idx_fm:int(nf / 2) + idx_fm + 1:, :] = 0

        fmin = fm
        fmax = f[-1]
    else:  # 'bandpass'
        fmin = float(filter_range.split('-')[0])
        fmax = float(filter_range.split('-')[1])

        mask1 = np.zeros((nf, nk), np.uint8)
        mask2 = np.zeros((nf, nk), np.uint8)
        scale_k = 1  # 3 / 4  # 0.25/5  # 缩放k区间比例
        if scale_k < 1:
            idx_fm = int(nf / fs * plot_derate * fmin)
            mask1[int(nf / 2) - idx_fm:int(nf / 2) + idx_fm + 1:,
            int(nk * (1 - scale_k) / 2):-int(nk * (1 - scale_k) / 2) + 1] = 1
            idx_fm = int(nf / fs * plot_derate * fmax)
            mask2[int(nf / 2) - idx_fm:int(nf / 2) + idx_fm + 1:,
            int(nk * (1 - scale_k) / 2):-int(nk * (1 - scale_k) / 2) + 1] = 1
        elif scale_k == 1:
            idx_fm = int(nf / fs * plot_derate * fmin)
            mask1[int(nf / 2) - idx_fm:int(nf / 2) + idx_fm + 1:, :] = 1
            idx_fm = int(nf / fs * plot_derate * fmax)
            mask2[int(nf / 2) - idx_fm:int(nf / 2) + idx_fm + 1:, :] = 1
        mask = mask2-mask1

    mFK_S = FK_S*mask
    # print(np.where(mFK_S!=0))
    filterd_img, dt, dx = Dt_by_2DiFFT(mFK_S, f[-1], k[-1])
    fig = plt.figure(dpi=300)
    plt.suptitle(file_path)
    plt.subplot(121)
    vmin2 = -3
    vmax2 = 3

    plot_data = stats.zscore(np.real(filterd_img), axis=None)

    # extent = [first_time_sample / fs, first_time_sample / fs+(plot_data[::, 0].shape[0] - 1) * plot_derate / fs,
    #           depth[first_channel], depth[last_channel]]
    img1 = plt.imshow(plot_data.T, cmap='bwr', aspect='auto', extent=extent, vmin=vmin2, vmax=vmax2,
                      origin='lower')
    plt.xlabel('Time (s)', fontsize=12)
    plt.ylabel('Distance (m)', fontsize=12)
    plt.xlim([40, 40.5])
    # plt.colorbar(img1)

    plt.plot([], color='w', linestyle="", label=f'{fmin}-{fmax} Hz')
    plt.legend()

    plt.subplot(122)
    plot_data = stats.zscore(np.abs(mFK_S), axis=None)
    img4 = plt.imshow(plot_data.T, cmap='bwr', aspect='auto', vmin=-vmax2, vmax=vmax2, extent=extent2,
                      origin='lower')
    plt.xlabel('Frequency (Hz)', fontsize=12)
    plt.ylabel('Wave number (1/m)', fontsize=12)
    plt.xlim([max(0, fmin-5),min(fmax+5,f[-1])])
    plt.colorbar(img4)
    plt.tight_layout()

    plt.savefig(f'filterd-{output_name}_{output_number}.png', format='png', bbox_inches='tight', pad_inches=0.1)
    plt.clf()
    # plt.show(block=False)

