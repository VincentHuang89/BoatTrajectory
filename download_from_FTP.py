import ftplib
import fnmatch
import os
from tqdm import tqdm



def is_ftp_connected(ftp):
    try:
        # 执行一个简单的FTP命令来检测连接状态
        ftp.voidcmd('NOOP')  # NOOP命令不执行任何操作，但可以用于检测连接状态
        return True  # 连接仍然有效

    except ftplib.all_errors:
        return False  # 连接已断开或不可用


def write_callback(file, data, pbar):
    file.write(data)
    data_size=len(data)
    pbar.update(data_size)


def download_file1(ftp, file, local_file):
    try:
        ftp.voidcmd('TYPE I')  # 设置传输模式为二进制模式
        total_size = ftp.size(file)
        if os.path.exists(local_file) and os.path.getsize(local_file) == total_size:
            str_info = print("file already exists.")
            return str_info
        if os.path.exists(local_file) and os.path.getsize(local_file) < total_size:  # 断点续传
            downloaded_size = os.path.getsize(local_file)
            with open(local_file, 'ab') as f:
                with tqdm(total=total_size, initial=downloaded_size, unit='B', unit_scale=True, desc=file) as pbar:
                    ftp.retrbinary('RETR ' + file, lambda data: write_callback(f, data, pbar), rest=downloaded_size)
            return 'Download completed.'
        else:
            with open(local_file, 'wb') as f:  # 完整下载
                with tqdm(total=total_size, unit='B', unit_scale=True, desc=file) as pbar:
                    ftp.retrbinary('RETR ' + file, lambda data: write_callback(f, data, pbar))
            return 'Download completed.'
    except ftplib.error_proto as e:
        # 捕获BrokenPipeError异常并进行处理
        if isinstance(e, BrokenPipeError):
            print("连接中断，尝试重新建立连接并重新下载文件...")
            ftp.quit()
            # 建立新的FTP连接
            ftp = create_FTP('172.18.218.194', 'hzd', 'optical503', '/NASShare4/sanjiao_guishan0904_1117/')  # 创建ftp对象并登录选择指定文件夹


            # 重新下载文件
            download_file(ftp, file, local_file)

        else:
            # 其他ftplib.error_proto异常的处理
            print("发生其他ftplib.error_proto异常:", e)

    except Exception as e:
        # 捕获其他异常并进行处理
        print("发生其他异常:", e)


def download_file(ftp, file, local_file):
    # ftp.set_pasv(True)  # 设置传输模式为被动模式(PASV模式)
    ftp.voidcmd('TYPE I')  # 设置传输模式为二进制模式
    total_size = ftp.size(file)
    if os.path.exists(local_file) and os.path.getsize(local_file) == total_size:
        str_info = print("file already exists.")
        return str_info
    if os.path.exists(local_file) and os.path.getsize(local_file) < total_size:  # 断点续传
        downloaded_size = os.path.getsize(local_file)
        with open(local_file, 'ab') as f:
            with tqdm(total=total_size, initial=downloaded_size, unit='B', unit_scale=True, desc=file) as pbar:
                ftp.retrbinary('RETR ' + file, lambda data: write_callback(f, data, pbar), rest=downloaded_size)
        return 'Download completed.'
    else:
        with open(local_file, 'wb') as f:  # 完整下载
            with tqdm(total=total_size, unit='B', unit_scale=True, desc=file) as pbar:
                ftp.retrbinary('RETR ' + file, lambda data: write_callback(f, data, pbar))
        return 'Download completed.'



# def download_file(ftp, file, local_file, MAX_RETRIES = 3,RETRY_DELAY = 20,CHUNK_SIZE = 100*1024*1024):
#     # MAX_RETRIES: 最大重传次数
#     # RETRY_DELAY:重试之间的延迟时间（秒）
#     # CHUNK_SIZE = 1024 * 1024  # 每个块的大小（字节）
#
#     total_size = ftp.size(file)
#     if os.path.exists(local_file) and os.path.getsize(local_file) == total_size:
#         str_info = print("file already exists.")
#         return str_info
#     if os.path.exists(local_file) and os.path.getsize(local_file) < total_size:  # 断点续传
#         downloaded_size = os.path.getsize(local_file)
#         with open(local_file, 'ab') as f:
#             with tqdm(total=total_size, initial=downloaded_size, unit='B', unit_scale=True, desc=file) as pbar:
#                 while downloaded_size < total_size:
#                     try:
#                         # 下载文件块
#                         ftp.retrbinary('RETR ' + file, lambda data: write_callback(f, data, pbar), rest=downloaded_size)
#                         downloaded_size = os.path.getsize(local_file)
#
#                         # 重置重试次数
#                         retries = 0
#
#                     except Exception as e:
#                         downloaded_size = os.path.getsize(local_file)
#                         print(f'下载失败: {str(e)}')
#                         retries += 1
#                         if retries >= MAX_RETRIES:
#                             raise
#
#                         # 等待一段时间后重试
#                         print(f'{RETRY_DELAY}秒后重试...')
#                         time.sleep(RETRY_DELAY)
#         return 'Download completed.'
#     else:
#         with open(local_file, 'wb') as f:  # 完整下载
#             with tqdm(total=total_size, unit='B', unit_scale=True, desc=file) as pbar:
#                 while True:
#                     try:
#                         # 下载文件块
#                         ftp.retrbinary('RETR ' + file, lambda data: write_callback(f, data, pbar))
#                         # 重置重试次数
#                         retries = 0
#
#                         break
#
#                     except Exception as e:
#                         print(f'下载失败: {str(e)}')
#                         retries += 1
#                         if retries >= MAX_RETRIES:
#                             raise
#
#                         # 等待一段时间后重试
#                         print(f'{RETRY_DELAY}秒后重试...')
#                         time.sleep(RETRY_DELAY)


def create_FTP(ip,user,password, dir, MAX_TIMEOUT=300):
    # 创建FTP对象
    ftp = ftplib.FTP()

    ftp.connect(ip, 21)
    # 登录FTP服务器
    ftp.login(user, password)

    # 设置传输模式为二进制模式
    ftp.sendcmd('TYPE I')

    # 切换到需要下载的目录
    ftp.cwd(dir)
    ftp.max_cons = 1  # 设置最大连接数为 2
    # ftp.set_pasv(True)  # 设置传输模式为被动模式(PASV模式)
    ftp.set_pasv(False)  # 设置传输模式为主动模式
    ftp.timeout =MAX_TIMEOUT
    # 设置数据连接超时
    # ftp.sock.settimeout(MAX_TIMEOUT)   # 最大等待时间（秒）
    return ftp


if __name__ == "__main__":
    # ftp = create_FTP('172.16.240.224', 'admin', 'optical503', '/0722/')  # 创建ftp对象并登录选择指定文件夹
    FTP_dir = '/NASShare4/sanjiao_guishan0904_1117/'
    #FTP_dir = ''

    ftp = create_FTP('172.18.218.194', 'hzd', 'optical503', FTP_dir)  # 创建ftp对象并登录选择指定文件夹
    #ftp = create_FTP('172.16.240.224', 'hzd', 'optical503', FTP_dir)  # 创建ftp对象并登录选择指定文件夹

    # 获取目录下的文件列表
    file_list = ftp.nlst()
    #print(file_list)
    # 指定保存下载文件的文件夹路径
    save_dir = r'/home/huangwj/DAS/BoatTrajectory/DataforAIS'
   

    '''遍历文件列表并下载符合条件的文件'''
    try:
        for file in file_list:
            # 删除数据
            # if fnmatch.fnmatch(file, '*20220722_06*.tdms'):   # 通配符匹配文件名
            #     local_file = os.path.join(save_dir, file)
            #     if os.path.exists(local_file):
            #         os.remove(local_file)
            #         print('delete:', file)

            # 下载数据
            if fnmatch.fnmatch(file, '*20211119_0110*.tdms'):   # 通配符匹配文件名
                print('Downloading:', file)
                local_file = os.path.join(save_dir, file)
                download_file(ftp, file, local_file)
    except:
        print('下载被中断！')
    finally:
        ftp.quit()  # 关闭FTP连接
