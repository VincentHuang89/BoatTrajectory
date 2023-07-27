'''
接口由协议、IP、端口、路径及参数组成，具体形式如下：

http://xxx.xxx.xxx.xxx:PORT/srvs?cmd=PARACMD&param=PARATEXT

岸基AIS：srvs  IALA：iala    卫星AIS：sate

上图是数据服务的通用接口形式，其中主要变化来自于参数部分：

cmd参数

请求命令码，不同的业务使用不同的请求命令码，具体参照数据服务的各自要求。

param参数

请求参数为JSON格式，具体参照数据API的各自要求。请求参数需使用base64进行编码。

注意：如果实时查询类服务的结果中某字段为null，则返回的json串中不会含有该字段。

'''


import requests
 
response = requests.get(url='http://httpbin.org/get',verify=False)
print(response.text)