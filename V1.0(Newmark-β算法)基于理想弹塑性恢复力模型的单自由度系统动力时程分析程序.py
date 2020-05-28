#Version 1.0 created by 罗星宇, 土木工学, 西南科技大学
#基于理想弹塑性恢复力模式和Newmak-β算法的单自由度系统动力时程分析
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

def NonlinReForce(fsi_1,Δu,F_srd0,k):#定义理想弹塑性恢复力-位移规则。fs(i-1) 为前一步的恢复力，Δu为位移差，k为初始刚度。
    if abs(fsi_1)<=F_srd0: #F_srd0 为初始屈服力(kN)。
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

################################      定义初始数据    ##################################
k=875.5#刚度系数
m=17.5#质量
c=35#粘滞阻尼系数
u0=0#初始位移
v0=0#初始速度
T=6#反应计算截至时刻
F_srd0=26.7#初始屈服力大小(kN)
γ=1/2
β=1/4
Δt=0.001#时间步长
N_Decimal=len(str(Δt).split(".")[1])#t保留小数位数与时间步长是对应的，二者应该同时修改。

if abs(u0*k)<F_srd0:#恢复力初值和初始位移有关。
    fs=u0*k
elif u0>0:
    fs=F_srd0
else:
    fs=-F_srd0
fs_ls=[]#用于存储恢复力之列表
fs_ls.append(fs)

u_ls=[]#定义位移列表，存储位移反应。
u_ls.append(u0)#位移初始条件u,0。
v_ls=[]#定义速度列表，存储速度反应。
v_ls.append(v0)#初始速度v0
a_ls=[]#定义加速度列表，存储加速度反应。
a_ls.append((P(0)-c*v0-fs)/m)#初始加速度a0
t=0#时程反应计时器
t_ls=[]#存储时程反应的时间散点
t_ls.append(0)
i=1#位移反应离散点个数。


################################      非线性时程分析   ###############################
t1=time.time()#非线性时程分析开始时刻t1
C1=m*pow(Δt,2)*(1-2*β)+c*pow(Δt,3)*(γ-2*β)
C2=2*β*pow(Δt,2)
C3=2*(m+c*γ*Δt)
C4=2*Δt*(m+c*Δt*(γ-β))
C5=2*(m+c*γ*Δt)
utemp=0#位移反应临时变量，最后存储到u列表中。


while(round(t,N_Decimal)<T):#当t=Ts时停止计算。完成本时间点的反应计算，再增加Δt，并对计数位移反应的个数i。
    t=t+Δt
    t_ls.append(t)
    utemp=(C1* a_ls[i-1]+C2*(P(t)-fs)+C3*u_ls[i-1]+C4*v_ls[i-1])/C5#当前的位移是用的前一点的恢复力计算得到的。这和线性的计算法有所区别。
    u_ls.append(utemp)
    fs=NonlinReForce(fs_ls[i-1],utemp-u_ls[i-1],F_srd0,k)#由于当前位移未知，所以用前一点位移的恢复力代替当前恢复力。
    fs_ls.append(fs)
    v_ls.append(γ*(u_ls[i]-u_ls[i-1])/(β*Δt)+(1-γ/β)*v_ls[i-1]+(1-γ/(2*β))*a_ls[i-1]*Δt)
    a_ls.append((u_ls[i]-u_ls[i-1])/(β*pow(Δt,2))-v_ls[i-1]/(β*Δt)-(1/(2*β)-1)*a_ls[i-1])
    i=i+1
t2=time.time()#非线性时程分析完成时刻t2

################################      绘制时程反应图    ###############################
t1_plot=time.time()#绘图过程初时刻t1
matplotlib.rcParams['figure.figsize'] = [15, 9.27] #画布的大小设置为15*9.27
matplotlib.rcParams['figure.subplot.left'] = 0.1#画布left间隙=10%
matplotlib.rcParams['figure.subplot.bottom'] = 0.1#画布bottom间隙=10%
matplotlib.rcParams['figure.subplot.right'] = 0.93#画布right间隙=7%
matplotlib.rcParams['figure.subplot.top'] = 0.93#画布top间隙=7%
plt.subplots_adjust(wspace =0.3, hspace =0.3)#调整子图间的宽度间隙，高度间隙。

plt.subplot(2,2,1)#2X2子图空间的第一个位置
plt.plot(t_ls,u_ls, 'black', label = 'the Newmark-β algorithm')#位移-时间曲线
plt.xlabel('时间(s)', fontdict={'family' : 'MicroSoft Yahei', 'size'   : 12})
plt.ylabel('位移(m)', fontdict={'family' : 'MicroSoft Yahei', 'size'   : 12})
plt.title("位移-时间曲线",fontdict={'family' : 'MicroSoft Yahei', 'size'   : 18,'weight': 'bold'})#图题
plt.legend(loc=1)#显示图例.loc=1的意思是将图例放在右上角
plt.yticks(fontproperties = 'Times New Roman', size = 10)
plt.xticks(fontproperties = 'Times New Roman', size = 10)

plt.subplot(2,2,2)#2X2子图空间的第二个位置
plt.plot(t_ls,v_ls, 'blue', label = 'the Newmark-β algorithm')#速度-时间曲线
plt.xlabel('时间(s)', fontdict={'family' : 'MicroSoft Yahei', 'size'   : 12})
plt.ylabel('速度(m/s)', fontdict={'family' : 'MicroSoft Yahei', 'size'   : 12})
plt.title("速度-时间曲线",fontdict={'family' : 'MicroSoft Yahei', 'size'   : 18,'weight': 'bold'})#图题
plt.legend(loc=1)#显示图例.loc=1的意思是将图例放在右上角
plt.yticks(fontproperties = 'Times New Roman', size = 10)
plt.xticks(fontproperties = 'Times New Roman', size = 10)

plt.subplot(2,2,3)#2X2子图空间的第三个位置
plt.plot(t_ls,a_ls, 'red', label = 'the Newmark-β algorithm')#加速度-时间曲线
plt.xlabel('时间(s)', fontdict={'family' : 'MicroSoft Yahei', 'size'   : 12})
plt.ylabel('加速度(m.s-2)', fontdict={'family' : 'MicroSoft Yahei', 'size'   : 12})
plt.title("加速度-时间曲线",fontdict={'family' : 'MicroSoft Yahei', 'size'   : 18,'weight': 'bold'})#图题
plt.legend(loc=1)#显示图例.loc=1的意思是将图例放在右上角
plt.yticks(fontproperties = 'Times New Roman', size = 10)
plt.xticks(fontproperties = 'Times New Roman', size = 10)


plt.subplot(2,2,4)#2X2子图空间的第四个位置
plt.plot(u_ls,fs_ls, 'purple', label = 'the Newmark-β algorithm')#力-位移曲线
plt.xlabel('位移(m)', fontdict={'family' : 'MicroSoft Yahei', 'size'   : 12})
plt.ylabel('恢复力(kN)', fontdict={'family' : 'MicroSoft Yahei', 'size'   : 12})
plt.title("恢复力-位移曲线",fontdict={'family' : 'MicroSoft Yahei', 'size'   : 18,'weight': 'bold'})#图题
plt.legend(loc=1)#显示图例.loc=1的意思是将图例放在右上角
plt.yticks(fontproperties = 'Times New Roman', size = 10)
plt.xticks(fontproperties = 'Times New Roman', size = 10)


t2_plot=time.time()#绘图过程末时刻t2
mngr = plt.get_current_fig_manager()  # 获取当前 画布 manager
mngr.window.wm_geometry("+205+20")  # 调整画布在屏幕上弹出的位置
plt.show()#显示图


################################      保存时程反应数据  ###############################
t_header="时间(s)"
u_header="位移(m)"
file_name_u_t="位移-时间数据.txt"
Output(t_ls,u_ls,t_header,u_header,file_name_u_t)#位移-时间数据
print('位移-时间数据已保存到"%s"！' %file_name_u_t) 

t_header="时间(s)"
v_header="速度(m/s)"
file_name_v_t="速度-时间数据.txt"
Output(t_ls,v_ls,t_header,v_header,file_name_v_t)#速度-时间数据
print('速度-时间数据已保存到"%s"！' %file_name_v_t) 

t_header="时间(s)"
a_header="加速度(m.s-2)"
file_name_a_t="加速度-时间数据.txt"
Output(t_ls,a_ls,t_header,a_header,file_name_a_t)#加速度-时间数据
print('加速度-时间数据已保存到"%s"！' %file_name_a_t) 


u_header="位移(m)"
fs_header="恢复力(kN)"
file_name_fs_u="恢复力-位移数据.txt"
Output(u_ls,fs_ls,u_header,fs_header,file_name_fs_u)#恢复力-位移数据
print('恢复力-位移数据已保存到"%s"！' %file_name_fs_u) 

print("时程分析耗时：{:.3f}s".format(t2-t1))#输出非线性分析耗时
print("绘图耗时：{:.3f}s".format(t2_plot-t1_plot))#输出绘图耗时
input()#等待输入，功能相当于暂停。防止窗口关闭，来不及查看信息。
################################      保存时程反应数据  ###############################

