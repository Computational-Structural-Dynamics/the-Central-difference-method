#Version 1.2 coded by Lo Hsingyu, 土木工学, 青义科技大学 
#单自由度系统，理想弹塑性动力时程分析(中心差分算法)
#单位：力kN，长度m，时间s，质量t，阻尼系数kNs/m，刚度kN/m
import numpy as np#用于生成等间距元素
import time#用于模块运行耗时计算
import matplotlib#用于修改画布的大小，及子图相对画布边缘的距离
from matplotlib import pyplot as plt#用于绘图


def P(t):#定义外部激励函数，0到0.4s，力从0变到40kN, 0.4s到0.8s,力从40kN变到0，呈现三角形变化。
     #return 40*np.sin(10*t)#testing exercite
    if 0<=t<=0.4:
        return 100*t
    elif 0.4<=t<=0.8:
        return 80-100*t
    else:
        return 0

def NonlinReForce(fsi_1,Δu,fs_srd,k):#理想弹塑性恢复力规则。fs(i-1) 前一点的恢复力，位移差Δu，k刚度
    if -fs_srd<=fsi_1<=fs_srd: #fs_srd 屈服力大小(kN)。
        if -fs_srd<=fsi_1+k*Δu<=fs_srd:
            return fsi_1+k*Δu
        elif 0<=fsi_1:
            return fs_srd
        else:
            return -fs_srd       
    elif 0<=fsi_1:
        return fs_srd
    else:
        return -fs_srd     

def Output(ls1,ls2,ls1_label,ls2_label,file_name):#定义将时程反应保存到文本的函数。按照ls1第一列，ls2第一列的顺序，将列表中的元素保存到文本。ls1_label和ls2_label分别是其表头，如t(s)和u(m)。
    file_write=open(file_name,"w+")
    file_write.write(str(ls1_label)+"   ")
    file_write.write(str(ls2_label)+"\n")
    for i in range(len(ls1)):
        ls1_temp=ls1[i]
        file_write.write(str(ls1_temp)+",")
        ls2_temp=ls2[i]
        file_write.write(str(ls2_temp)+"\n")
    file_write.close()  

################################      定义初始数据    ##################################
k=875.5#刚度系数
m=17.5#质量
c=35#粘滞阻尼系数
T=1.2#反应计算截至时刻
fs_srd=26.7#屈服力大小
Δt=0.01#时间步长
N_Decimal=2#t保留小数位数与时间步长是对应的，二者必须同时修改。
u_ls=[]#定义位移列表，存储位移反应。
u_ls.append(0)#位移初始条件u,-1。
u_ls.append(0)#位移初始条件u,0。
v_ls=[]#定义速度列表，存储速度反应。
v_ls.append(0)
a_ls=[]#定义加速度列表，存储加速度反应。
a_ls.append(0)
t=0
i=2#位移反应离散点个数，因为已经有了2个初始条件u[0]=u[1]=0，所以计数从2开始。反应计算完成后需要删除第一个元素。


################################      非线性时程分析   ###############################
t1=time.time()#中心差分算法，非线性时程分析初时刻t1
utemp=0#临时变量，将每一次计算出的位移反应。最后存储到u列表中
fs=0#恢复力,初值为0。
fs_ls=[]#用于保存恢复力的列表
fs_ls.append(0)
while(round(t,N_Decimal)<T):#当t=Ts时停止计算。完成本时间点的反应计算，再增加Δt，并对计数位移反应的个数i。
    i=i+1#先对位移进行计数再计算位移。例如先将i-->5，再计算第5个反应。
    utemp=(P(t)-fs+2*m/pow(Δt,2)*u_ls[i-2]-(m/pow(Δt,2)-c/(2*Δt))*u_ls[i-3])/(m/pow(Δt,2)+c/(2*Δt))#第i个位移反应
    t=t+Δt
    u_ls.append(utemp)
    v_ls.append((u_ls[i-1]-u_ls[i-3])/(2*Δt))
    a_ls.append((u_ls[i-1]-2*u_ls[i-2]+u_ls[i-3])/pow(Δt,2))
    fs=NonlinReForce(fs,u_ls[i-1]-u_ls[i-2],fs_srd,k)#第i个恢复力
    fs_ls.append(fs)
del u_ls[0]#删除第一个元素，即u,-1,因为u,-1不是时程反应中的数据，仅仅是2个起步位移之一。
t_ls=np.arange(0,T+Δt,Δt)#数值解的时刻列表，0到T，间隔Δt形成离散点。
t2=time.time()#中心差分算法，非线性时程分析末时刻t2

    
################################      绘制时程反应图    ###############################
t1_plot=time.time()#绘图过程初时刻t1


matplotlib.rcParams['figure.figsize'] = [15, 9.27] #画布的大小设置为15*9.27
matplotlib.rcParams['figure.subplot.left'] = 0.1#画布left间隙=10%
matplotlib.rcParams['figure.subplot.bottom'] = 0.1#画布bottom间隙=10%
matplotlib.rcParams['figure.subplot.right'] = 0.93#画布right间隙=7%
matplotlib.rcParams['figure.subplot.top'] = 0.93#画布top间隙=7%
plt.subplots_adjust(wspace =0.3, hspace =0.3)#调整子图间的宽度间隙，高度间隙。


fig=plt.subplot(2,2,1)#2X2子图空间的第一个位置
plt.plot(t_ls,u_ls, 'black', label = 'the Central difference algorithm')#位移时间曲线
plt.xlabel('时间(s)', fontdict={'family' : 'MicroSoft Yahei', 'size'   : 12})
plt.ylabel('位移(m)', fontdict={'family' : 'MicroSoft Yahei', 'size'   : 12})
plt.title("位移-时间曲线",fontdict={'family' : 'MicroSoft Yahei', 'size'   : 18,'weight': 'bold'})#图题
plt.legend()#显示图例
plt.yticks(fontproperties = 'Times New Roman', size = 10)
plt.xticks(fontproperties = 'Times New Roman', size = 10)

plt.subplot(2,2,2)#2X2子图空间的第二个位置
plt.plot(t_ls,v_ls, 'blue', label = 'the Central difference algorithm')#速度时间曲线
plt.xlabel('时间(s)', fontdict={'family' : 'MicroSoft Yahei', 'size'   : 12})
plt.ylabel('速度(m/s)', fontdict={'family' : 'MicroSoft Yahei', 'size'   : 12})
plt.title("速度-时间曲线",fontdict={'family' : 'MicroSoft Yahei', 'size'   : 18,'weight': 'bold'})#图题
plt.legend()#显示图例
plt.yticks(fontproperties = 'Times New Roman', size = 10)
plt.xticks(fontproperties = 'Times New Roman', size = 10)

plt.subplot(2,2,3)#2X2子图空间的第三个位置
plt.plot(t_ls,a_ls, 'red', label = 'the Central difference algorithm')#加速度时间曲线
plt.xlabel('时间(s)', fontdict={'family' : 'MicroSoft Yahei', 'size'   : 12})
plt.ylabel('加速度(m.s-2)', fontdict={'family' : 'MicroSoft Yahei', 'size'   : 12})
plt.title("加速度-时间曲线",fontdict={'family' : 'MicroSoft Yahei', 'size'   : 18,'weight': 'bold'})#图题
plt.legend()#显示图例
plt.yticks(fontproperties = 'Times New Roman', size = 10)
plt.xticks(fontproperties = 'Times New Roman', size = 10)

plt.subplot(2,2,4)#2X2子图空间的第四个位置
plt.plot(u_ls,fs_ls, 'purple', label = 'the Single degree system')#力位移曲线
plt.xlabel('位移(m)', fontdict={'family' : 'MicroSoft Yahei', 'size'   : 12})
plt.ylabel('恢复力(kN)', fontdict={'family' : 'MicroSoft Yahei', 'size'   : 12})
plt.title("恢复力-位移曲线",fontdict={'family' : 'MicroSoft Yahei', 'size'   : 18,'weight': 'bold'})#图题
plt.legend()#显示图例
plt.yticks(fontproperties = 'Times New Roman', size = 10)
plt.xticks(fontproperties = 'Times New Roman', size = 10)
t2_plot=time.time()#绘图过程末时刻t2

mngr = plt.get_current_fig_manager()  # 获取当前 画布 manager
mngr.window.wm_geometry("+205+20")  # 调整画布在屏幕上弹出的位置
plt.show()#显示图



################################      保存时程反应数据  ###############################
t_header="时间(s)"
u_header="位移(m)"
file_name_u_t="中心差分算法_位移-时间.txt"
Output(t_ls,u_ls,t_header,u_header,file_name_u_t)#位移-时间数据
print('位移-时间数据已保存到"%s"！' %file_name_u_t) 

t_header="时间(s)"
v_header="速度(m/s)"
file_name_v_t="中心差分算法_速度-时间.txt"
Output(t_ls,v_ls,t_header,v_header,file_name_v_t)#速度-时间数据
print('速度-时间数据已保存到"%s"！' %file_name_v_t) 

t_header="时间(s)"
a_header="加速度(m.s-2)"
file_name_a_t="中心差分算法_加速度-时间.txt"
Output(t_ls,a_ls,t_header,a_header,file_name_a_t)#加速度-时间数据
print('加速度-时间数据已保存到"%s"！' %file_name_a_t) 

u_header="位移(m)"
fs_header="恢复力(kN)"
file_name_fs_u="中心差分算法_恢复力-位移.txt"
Output(u_ls,fs_ls,u_header,fs_header,file_name_fs_u)#恢复力-位移数据
print('恢复力-位移数据已保存到"%s"！' %file_name_fs_u) 

print("非线性时程分析耗时：{:.3f}s".format(t2-t1))#输出非线性分析耗时
print("绘图耗时：{:.3f}s".format(t2_plot-t1_plot))#输出绘图耗时
input()#等待输入，功能相当于暂停。
################################      保存时程反应数据  ###############################
