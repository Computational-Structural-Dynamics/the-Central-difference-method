#Version 1.4 created by 罗星宇, 土木工学, 西南科技大学
#基于理想弹塑性恢复力模式和中心差分算法的单自由度系统动力时程分析
#单位：力kN，长度m，时间s，质量t，粘滞阻尼系数kNs/m，刚度系数kN/m
import time#用于模块运行耗时计算
import matplotlib#用于修改画布的大小和子图相对画布边缘的距离
import numpy as np
from matplotlib import pyplot as plt#用于绘图

def P(t):#定义外部激励，0到0.4s，力从0变到40kN, 0.4s到0.8s,力从40kN变到0，呈现三角形变化。
    #return 0

    if 0<=t<=0.4:
        return 100*t
    elif 0.4<t<=0.8:
        return 80-100*t
    else:
        return 0

    
def ReForce(fsi_1,Δu,F_srd0,k):#定义理想弹塑性恢复力-位移规则。fs(i-1) 为前一步的恢复力，Δu为位移差，k为初始刚度。
    if abs(fsi_1)<=F_srd0: #F_srd0为屈服力大小(kN)。
        if abs(fsi_1+k*Δu)<=F_srd0:
            return fsi_1+k*Δu
        elif 0<fsi_1+k*Δu:
            return F_srd0
        else:
            return -F_srd0       
    elif fsi_1>0:
        if Δu>0:
            return F_srd0
        else:
            return fsi_1+k*Δu
    else:
        if Δu<0:
            return -F_srd0
        else:
            return fsi_1+k*Δu      

def Output(ls1,ls2,ls1_header,ls2_header,file_name):#定义保存时程反应。按照ls1第一列，ls2第一列的顺序，将列表中的元素保存到文本。ls1_header和ls2_header分别是其表头，如时间(s)和位移(m)。
    file_write=open(file_name,"w+")
    file_write.write(str(ls1_header)+"   ")
    file_write.write(str(ls2_header)+"\n")
    for i in range(len(ls1)):
        file_write.write(str(ls1[i])+",")
        file_write.write(str(ls2[i])+"\n")
    file_write.close()

def Figset(plotsize):
    matplotlib.rcParams['figure.figsize'] =plotsize #设置画布大小
    matplotlib.rcParams['figure.subplot.left'] = 0.1#画布left间隙=10%
    matplotlib.rcParams['figure.subplot.bottom'] = 0.1#画布bottom间隙=10%
    matplotlib.rcParams['figure.subplot.right'] = 0.93#画布right间隙=7%
    matplotlib.rcParams['figure.subplot.top'] = 0.93#画布top间隙=7%
    plt.subplots_adjust(wspace =0.3, hspace =0.3)#调整子图间的宽度间隙，高度间隙。
    mngr = plt.get_current_fig_manager()  # 获取当前 画布 manager
    mngr.window.wm_geometry("+205+20")  # 调整画布在屏幕上弹出的位置
class Datamngr:
    def __init__(self,dmtn1,dmtn2):
        self.dmtn1=dmtn1      
        self.dmtn2=dmtn2        
    def plot(self,place,ls1,ls2,color,xlabe,ylabel,title):
        plt.subplot(self.dmtn1,self.dmtn2,place)#2X2子图空间的第二个位置
        plt.plot(ls1,ls2, color, label = 'the Central difference algorithm')#速度-时间曲线
        plt.xlabel(xlabe, fontdict={'family' : 'MicroSoft Yahei', 'size'   : 12})
        plt.ylabel(ylabel, fontdict={'family' : 'MicroSoft Yahei', 'size'   : 12})
        plt.title(title,fontdict={'family' : 'MicroSoft Yahei', 'size'   : 18,'weight': 'bold'})#图题
        plt.legend(loc=1)#显示图例.loc=1的意思是将图例放在右上角
        plt.yticks(fontproperties = 'Times New Roman', size = 10)
        plt.xticks(fontproperties = 'Times New Roman', size = 10)
    def Savedata(self,ls1,ls2,ls1_header,ls2_header,file_nm):
        Output(ls1,ls2,ls1_header, ls2_header,file_nm)#速度-时间数据
        print('"%s"已保存！' %file_nm) 
################################      定义初始数据    ##################################
k=875.5#刚度系数
m=17.5#质量
c=35#粘滞阻尼系数
u0=0#初始位移
v0=0#初始速度
T=6#反应计算截至时刻
F_srd0=26.7#屈服力大小
Δt=0.001#时间步长
N_Decimal=3#t保留小数位数与时间步长是对应的，二者应该同时修改。

fs=ReForce(0,u0,F_srd0,k)#恢复力初值和初始位移及初始刚度系数有关。
fs_ls=[]#用于存储恢复力之列表
fs_ls.append(fs)

u_ls=[]#定义位移列表，存储位移反应。
u_ls.append(u0-Δt*v0+pow(Δt,2)/2*(P(0)-c*v0-k*u0)/m)#位移初始条件u,-1。
u_ls.append(u0)#位移初始条件u,0。
v_ls=[]#定义速度列表，存储速度反应。
v_ls.append(v0)#初始速度v0
a_ls=[]#定义加速度列表，存储加速度反应。
a_ls.append((P(0)-c*v0-fs)/m)#初始加速度a0
t=0#时程反应计时器
t_ls=[]#存储时程反应的时间散点
t_ls.append(0)
i=2#位移反应离散点个数，因为已经有了2个初始条件u[0]=u[1]=0，所以计数从2开始。完成时程反应计算后，需要删除第一个元素。


################################      非线性时程分析   ###############################
t1=time.time()#非线性时程分析开始时刻t1
utemp=0#位移反应临时变量，最后存储到u列表中。

while(round(t,N_Decimal)<T):#当t=Ts时停止计算。完成本时间点的反应计算，再增加Δt，并对计数位移反应的个数i。
    i=i+1#先对位移进行计数再计算位移。例如先将i-->5，再计算第5个反应。
    t=t+Δt#在0初始条件下(u0=0,v0=0), 为了不使等式右侧各项都等于零，所以用当前时刻的t代替前一步的t,因为Δt足够小，所以误差得到控制。
    #在这一点上与书上公式不同。书上的P(t)用的是前一步长的t, 这里用的是当前时刻的t。
    t_ls.append(round(t,N_Decimal))
    utemp=(P(t)-fs+2*m/pow(Δt,2)*u_ls[i-2]-(m/pow(Δt,2)-c/(2*Δt))*u_ls[i-3])/(m/pow(Δt,2)+c/(2*Δt))
    u_ls.append(utemp)
    v_ls.append((u_ls[i-1]-u_ls[i-3])/(2*Δt))
    a_ls.append((u_ls[i-1]-2*u_ls[i-2]+u_ls[i-3])/pow(Δt,2))
    fs=ReForce(fs,u_ls[i-1]-u_ls[i-2],F_srd0,k)#第i个恢复力
    fs_ls.append(fs)
del u_ls[0]#删除第一个元素，即u,-1,因为u,-1不是时程反应中的数据，仅仅是2个起步位移之一。
t2=time.time()#非线性时程分析完成时刻t2

    
################################      绘制时程反应图    ###############################
t1_plot=time.time()#绘图过程初时刻t1
dmtn1=2#子图行数
dmtn2=3#子图列数
plotsize=[15,9.27]
t_label='时间(s)'
u_label='位移(m)'
v_label='速度(m/s)'
a_label='加速度(m.s-2)'
fs_label='恢复力(kN)'

u_t_title="位移-时间曲线"
v_t_title="速度-时间曲线"
a_t_title="加速度-时间曲线"
fs_t_title="恢复力-时间曲线"
fs_u_title="恢复力-位移曲线"

u_t_color="black"
v_t_color="blue"
a_t_color="pink"
fs_t_color="red"
fs_u_color="grey"

u_t_place=1
v_t_place=2
a_t_place=3
fs_t_place=4
fs_u_place=5

Figset(plotsize)#设置画布大小准备绘图
u_t_Datamngr=Datamngr(dmtn1,dmtn2,)#创建u_t_数据操作器
u_t_Datamngr.plot(u_t_place,t_ls,u_ls,u_t_color,t_label,u_label,u_t_title)
v_t_Datamngr=Datamngr(dmtn1,dmtn2,)#创建v_t_数据操作器
v_t_Datamngr.plot(v_t_place,t_ls,v_ls,v_t_color,t_label,v_label,v_t_title)
a_t_Datamngr=Datamngr(dmtn1,dmtn2,)#创建a_t_数据操作器
a_t_Datamngr.plot(a_t_place,t_ls,a_ls,a_t_color,t_label,a_label,a_t_title)
fs_t_Datamngr=Datamngr(dmtn1,dmtn2,)#创建fs_t_数据操作器
fs_t_Datamngr.plot(fs_t_place,t_ls,fs_ls,fs_t_color,t_label,fs_label,fs_t_title)
fs_u_Datamngr=Datamngr(dmtn1,dmtn2,)#创建fs_u_数据操作器
fs_u_Datamngr.plot(fs_u_place,u_ls,fs_ls,fs_u_color,u_label,fs_label,fs_u_title)
t2_plot=time.time()#绘图过程末时刻t2
plt.show()


################################      保存时程反应数据  ###############################
file_name_u_t="位移-时间数据.txt"
file_name_v_t="速度-时间数据.txt"
file_name_a_t="加速度-时间数据.txt"
file_name_fs_t="恢复力-时间数据.txt"
file_name_fs_u="恢复力-位移数据.txt"

u_t_Datamngr.Savedata(t_ls,u_ls,t_label,u_label,file_name_u_t)#保存数据
v_t_Datamngr.Savedata(t_ls,v_ls,t_label,v_label,file_name_v_t)
a_t_Datamngr.Savedata(t_ls,a_ls,t_label,a_label,file_name_a_t)
fs_t_Datamngr.Savedata(t_ls,fs_ls,t_label,fs_label,file_name_fs_t)
fs_u_Datamngr.Savedata(u_ls,fs_ls,u_label,fs_label,file_name_fs_u)

print("时程分析耗时：{:.3f}s".format(t2-t1))#输出非线性分析耗时
print("绘图耗时：{:.3f}s".format(t2_plot-t1_plot))#输出绘图耗时
input()#等待输入，功能相当于暂停。防止窗口关闭，来不及查看信息。
################################      保存时程反应数据  ###############################

