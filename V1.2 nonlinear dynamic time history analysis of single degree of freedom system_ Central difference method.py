#Version 1.2 created by Lo Hsingyu, SWUST
#SDOF，force-displacement curve is based on the Ideal Elastoplasticity
#using the Central difference algorithm
#Units：force:kN，length:meter，time:second，mass:ton，viscous damping coefficient:kNs/m，stiffness coefficient:kN/m
import numpy as np#Used to generate list elements
import time#Used for calculation of module running time
import matplotlib#Used to modify the size of the canvas and the distance between the subgraph and the edge of the canvas
from matplotlib import pyplot as plt#drawing

def P(t):#define excitation
    if 0<=t<=0.4:
        return 100*t
    elif 0.4<=t<=0.8:
        return 80-100*t
    else:
        return 0

def NonlinReForce(fsi_1,Δu,fs_srd,k):#the Ideal Elastoplasticity。fs (i-1) is the resilience of the previous point. Δ u is the displacement difference. K is stiffness
    if -fs_srd<=fsi_1<=fs_srd: #fs_ srd is the yield force (KN).
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

def Output(ls1,ls2,ls1_label,ls2_label,file_name):#Defines a function that saves a time history response to text. Save the elements in the list to text in the order of LS1 first column and Ls2 first column. ls1_ Label and Ls2_ Label is its header, such as t (s) and U (m).
    file_write=open(file_name,"w+")
    file_write.write(str(ls1_label)+"   ")
    file_write.write(str(ls2_label)+"\n")
    for i in range(len(ls1)):
        ls1_temp=ls1[i]
        file_write.write(str(ls1_temp)+",")
        ls2_temp=ls2[i]
        file_write.write(str(ls2_temp)+"\n")
    file_write.close()  

################################      Define initial data    ##################################
k=875.5#Stiffness coefficient
m=17.5#mass
c=35#viscous damping coefficient
T=1.2#calculating duration
fs_srd=26.7#yield force
Δt=0.01#time step
N_Decimal=2#keep to how many decimal places. corresponding to Δt.
u_ls=[]#define displacement list and it will be used to store displacement response.
u_ls.append(0)#Initial condition of displacement.u,-1。
u_ls.append(0)#Initial condition of displacement.u,0。
v_ls=[]#define velocity list and it will be used to store velocity response.
v_ls.append(0)
a_ls=[]#Define acceleration list and it will be used to store acceleration response.
a_ls.append(0)
t=0
i=2#The displacement reflects the number of discrete points. Since there are already two initial conditions u [0] = u [1] = 0, the count starts from 2. The first element needs to be removed after reaction calculation.


################################      Nonlinear time history analysis   ###############################
t1=time.time()#Central difference algorithm. Initial time t1 of nonlinear time history analysis
utemp=0#A temporary variable will be calculated for each displacement response. Last stored inu list
fs=0#Resilience,its initial value is 0.
fs_ls=[]#Resilience list and it will be used to store resilience.
fs_ls.append(0)
while(round(t,N_Decimal)<T):#Stop calculation when t = ts. After the reaction calculation at this time point is completed, add Δ T, and count the number of displacement reactions I.
    i=i+1#First count the displacement and then calculate the displacement. For example, first i -- > 5, and then calculate the fifth reaction.
    utemp=(P(t)-fs+2*m/pow(Δt,2)*u_ls[i-2]-(m/pow(Δt,2)-c/(2*Δt))*u_ls[i-3])/(m/pow(Δt,2)+c/(2*Δt))#The ith displacement response
    t=t+Δt
    u_ls.append(utemp)
    v_ls.append((u_ls[i-1]-u_ls[i-3])/(2*Δt))
    a_ls.append((u_ls[i-1]-2*u_ls[i-2]+u_ls[i-3])/pow(Δt,2))
    fs=NonlinReForce(fs,u_ls[i-1]-u_ls[i-2],fs_srd,k)#the ith resilience
    fs_ls.append(fs)
del u_ls[0]#Delete the first element, u, - 1, because u, - 1 is not the data in the time history reaction, only one of the two start displacements.
t_ls=np.arange(0,T+Δt,Δt)#The time list of numerical solution, 0 to t, interval Δ t form discrete points.
t2=time.time()# t2 at the end of nonlinear time history analysis

    
################################      Drawing time history response fig    ###############################
t1_plot=time.time()#t1. start of drawing


matplotlib.rcParams['figure.figsize'] = [15, 9.27] #set canvas size  to 15 * 9.27
matplotlib.rcParams['figure.subplot.left'] = 0.1#Canvas left gap = 10%
matplotlib.rcParams['figure.subplot.bottom'] = 0.1#Canvas bottom gap = 10%
matplotlib.rcParams['figure.subplot.right'] = 0.93#Canvas right gap = 7%
matplotlib.rcParams['figure.subplot.top'] = 0.93#Canvas top gap = 7%
plt.subplots_adjust(wspace =0.3, hspace =0.3)#Adjust the width gap and height gap between the subplots.


fig=plt.subplot(2,2,1)#The first position of 2x2 subplots space
plt.plot(t_ls,u_ls, 'black', label = 'the Central difference algorithm')#the displacement-time curve
plt.xlabel('Time(s)', fontdict={'family' : 'MicroSoft Yahei', 'size'   : 12})
plt.ylabel('Displacement(m)', fontdict={'family' : 'MicroSoft Yahei', 'size'   : 12})
plt.title("Displacement-Time Curve",fontdict={'family' : 'MicroSoft Yahei', 'size'   : 18,'weight': 'bold'})#the title of single subplot
plt.legend()#legend of this subplot
plt.yticks(fontproperties = 'Times New Roman', size = 10)
plt.xticks(fontproperties = 'Times New Roman', size = 10)

plt.subplot(2,2,2)#The second position of 2x2 subplots space
plt.plot(t_ls,v_ls, 'blue', label = 'the Central difference algorithm')#the velocity-time curve
plt.xlabel('Time(s)', fontdict={'family' : 'MicroSoft Yahei', 'size'   : 12})
plt.ylabel('Velocity(m/s)', fontdict={'family' : 'MicroSoft Yahei', 'size'   : 12})
plt.title("Velocity-Time Curve",fontdict={'family' : 'MicroSoft Yahei', 'size'   : 18,'weight': 'bold'})#the title of single subplot
plt.legend()#legend of this subplot
plt.yticks(fontproperties = 'Times New Roman', size = 10)
plt.xticks(fontproperties = 'Times New Roman', size = 10)

plt.subplot(2,2,3)#The third position of 2x2 subplots space
plt.plot(t_ls,a_ls, 'red', label = 'the Central difference algorithm')#the acceleration-time curve
plt.xlabel('Time(s)', fontdict={'family' : 'MicroSoft Yahei', 'size'   : 12})
plt.ylabel('Acceleration(m.s-2)', fontdict={'family' : 'MicroSoft Yahei', 'size'   : 12})
plt.title("Acceleration-Time Curve",fontdict={'family' : 'MicroSoft Yahei', 'size'   : 18,'weight': 'bold'})#the title of single subplot
plt.legend()#legend of this subplot
plt.yticks(fontproperties = 'Times New Roman', size = 10)
plt.xticks(fontproperties = 'Times New Roman', size = 10)

plt.subplot(2,2,4)#The fourth position of 2x2 subplots space
plt.plot(u_ls,fs_ls, 'purple', label = 'the Single degree system')#the Resilience-displacement curve
plt.xlabel('Displacement(m)', fontdict={'family' : 'MicroSoft Yahei', 'size'   : 12})
plt.ylabel('Resilience(kN)', fontdict={'family' : 'MicroSoft Yahei', 'size'   : 12})
plt.title("Resilience-Displacement Curve",fontdict={'family' : 'MicroSoft Yahei', 'size'   : 18,'weight': 'bold'})#the title of single subplot
plt.legend()#legend of this subplot
plt.yticks(fontproperties = 'Times New Roman', size = 10)
plt.xticks(fontproperties = 'Times New Roman', size = 10)
t2_plot=time.time()#responses process ends ,t2

mngr = plt.get_current_fig_manager()  #get current canvas manager
mngr.window.wm_geometry("+205+20")  # Adjust where the canvas pops up on the screen
plt.show()#show the figure, including four subplots



################################      保存时程反应数据  ###############################
t_header="Time(s)"
u_header="Displacement(m)"
file_name_u_t="Displacement-Time Data.txt"
Output(t_ls,u_ls,t_header,u_header,file_name_u_t)#Displacement-Time Data
print('Displacement-Time Data has been saved to"%s"！' %file_name_u_t) 

t_header="Time(s)"
v_header="Velocity(m/s)"
file_name_v_t="Velocity-Time Data.txt"
Output(t_ls,v_ls,t_header,v_header,file_name_v_t)#Velocity-Time Data
print('Velocity-Time Data has been saved to"%s"！' %file_name_v_t) 

t_header="Time(s)"
a_header="Acceleration(m.s-2)"
file_name_a_t="Acceleration-Time Data.txt"
Output(t_ls,a_ls,t_header,a_header,file_name_a_t)#Acceleration-Time Data
print('Acceleration-Time Data has been saved to"%s"！' %file_name_a_t) 

u_header="Displacement(m)"
fs_header="Resilience(kN)"
file_name_fs_u="Displacement-Resilience Data.txt"
Output(u_ls,fs_ls,u_header,fs_header,file_name_fs_u)#Displacement-Resilience Data
print('Displacement-Resilience Data has been saved to"%s"！' %file_name_fs_u) 

print("Analysis process takes {:.3f}s".format(t2-t1))#Time consuming in nonlinear time history analysis
print("Drawing process takes {:.3f}s".format(t2_plot-t1_plot))#Time consuming in drawing
input()#Wait for input, the function is equivalent to pause.
################################      保存时程反应数据  ###############################
