import matplotlib.pyplot as plt
import pandas as pd

def dx(CO2_Air,CO2_Top):
    ret=[CO2_Air,CO2_Top]
    ret[0]= (MC_BlowAir() + MC_ExitAir() + MC_PadAir(CO2_Air)-MC_AirCan(CO2_Air)-MC_AirTop(CO2_Air,CO2_Top)-MC_TopOut(CO2_Top)-MC_AirOut(CO2_Air))/capAir
    ret[1]= (MC_AirTop(CO2_Air,CO2_Top)-MC_TopOut(CO2_Top))/capTop
    #print ("dx:")
    #print(ret)
    return ret
# RK-4 method
def rk4(CO2Air,CO2Top,h):
    
    # create list for return
    ret = [CO2Air,CO2Top]
    # create var to hold dx value 
    #Runge-kutta code
    k1=h*dx(ret[0],ret[1])
    k2=h*dx(ret[0]+k1[0]/2,ret[1]+k1[1]/2)
    k3=h*dx(ret[0]+k2[0]/2,ret[1]+k2[1]/2)
    k4=h*dx(ret[0]+k3[0]  ,ret[1]+k3[1]  )
    k=[0,0]
    i=0
    while i < 2:
        k[i]=(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6
        ret[i]+=k[i]
        i+=1
    #print("k= ")
    #print(k)
    return ret

def euler(CO2_Air,CO2_Top,h):
    
    ret= [CO2_Air,CO2_Top]
    #print("slope")
    # print(ret[0],ret[1])
    # print("Ret:",ret)
    slope = dx(ret[0],ret[1])
    # print("Slope: ",slope)
    i=0
    while i < 2:
        ret[i]+=slope[i]*h
        i+=1
    return ret
# we assume our result func have the format of math.exp(kx)

import math


    # 3           mg/J    [0,1]    N/A    m2
def MC_BlowAir():
    ans = n_HeatCO2 * U_Blow * P_Blow
    ans /= A_Flr

    return ans

    # 4          N/A       mg/s     m2


def MC_ExitAir():
    ans = U_ExtCO2 * o_ExtCO2
    ans /= A_Flr

    return ans

    # 5        m/s     mg/s     mg/s


def MC_PadAir(CO2Air):
    ans = f_Pad() * (CO2_Out - CO2Air)

    return ans

    # 5.1 [0,1]   m3/s    m2


def f_Pad():
    ans = U_Pad * o_Pad
    ans /= A_Flr

    return ans

    # 6         m/s     mg/s     mg/s


def MC_AirTop(CO2Air,CO2Top):
    ans = f_ThScr() * (CO2Air - CO2Top)

    return ans

    # 7        N/A || m(K^2/3)/s  K     K   m/s2  kg/m3     =      =


def f_ThScr():
    ans = U_ThScr * K_ThScr * (abs(T_Air - T_Top)) ** (2 / 3)
    ans += (1 - U_ThScr) * ((g * (1 - U_ThScr)) / (2 * p_MeanAir) * abs(p_Air - p_Top)) ** (1 / 2)

    return ans

    # 8      m   m   kg/m3|9.81|kg/m3


def o_crack(L, SO, p_Mean, p1, p2):
    ans = L * SO / p_Mean
    ans *= ((1 / 2) * p_Mean * SO * g * (p1 - p2)) ** (1 / 2)

    return ans

    # 9   khong su dung phuong trinh nay, day la ham cua nghien cuu cua [MG] va cong su de giai quyet navier stoke


def MC_AirOut(CO2Air):
    ans = (f_VentSide() + f_VentForced(o_VentForced)) * (CO2Air - CO2_Out)

    return ans

    # 10            N/A    m2                     m2      m2        m        K      K        K      N/A    m/s


def f_VentRoofSide(A_Roof):
    ans = C_d / A_Flr

    temp1 = (U_Roof ** 2 * U_Side ** 2 * A_Roof ** 2 * A_Side ** 2) / (
            U_Roof ** 2 * A_Roof ** 2 + U_Side ** 2 * A_Side ** 2)
    temp2 = (2 * g * h_SideRoof * (T_Air - T_Out)) / T_MeanAir
    temp3 = ((U_Roof * A_Roof + U_Side * A_Side) / 2) ** 2 * C_w * v_Wind ** 2

    ans *= (temp1 * temp2 + temp3) ** (1 / 2)

    return ans

    # 11        N/A


def n_InsScr():
    ans = S_InsScr * (2 - S_InsScr)

    return ans

    # 12        m/s       N/A


def f_leakage():
    if v_Wind < 0.25:
        ans = 0.25 * c_leakage
    else:
        ans = v_Wind * c_leakage

    return ans

    # 13
def f2_VentSide():
    ans = C_d*U_Side*A_Side*v_Wind/(2*A_Flr)*(C_w**(1/2))
    return ans

def f_VentSide():
    if n_Side >= n_Side_Thr:
        ans = n_InsScr() * f2_VentSide() + 0.5 * f_leakage()
    else:
        ans = n_InsScr() * (U_ThScr * f2_VentSide() + (1 - U_ThScr) * f_VentRoofSide(A_Roof) * n_Side) + 0.5 * f_leakage()

    return ans

    # 14


def f_VentForced(o_VentForced):
    ans = n_InsScr() * U_VentForced * o_VentForced
    ans /= A_Flr

    return ans

    # 15           m/s       mg/s     mg/s


def MC_TopOut(CO2Top):
    ans = f_VentRoof() * (CO2Top - CO2_Out)

    return ans

    # 16


def f_VentRoof():
    if n_Roof > n_Roof_Thr:
        ans = n_InsScr() * f2_VentRoof() + 0.5 * f_leakage()
    else:
        ans = n_InsScr() * (U_ThScr * f2_VentRoof() + (1 - U_ThScr) * f_VentRoofSide(A_Roof) + 0.5 * f_leakage())
    return ans

    # 17


def f2_VentRoof():
    ans = (C_d * U_Roof * A_Roof) / (2 * A_Flr)
    temp1 = (g * h_Roof * (T_Air - T_Out)) / (2 * T_MeanAir)
    temp2 = C_w * v_Wind ** 2

    ans =ans *( abs(temp1 + temp2) ** (1 / 2))

    return ans

    # 18     mg/mcmol     mg/m2        mcmol/m2.s

    #new cannabis section
def JPot():
    ans=JMax * (math.exp(Ej*(T_Can-T_25)/(R*T_Can*T_25))) * (1+math.exp((S*T_25-H)/(R*T_25)))/(1+math.exp((S*T_Can-H)/(R*T_Can)))
    return ans
def PAR_GHCan():
    ans = PAR_GH()*(1-p_Can)*(1-math.exp(-K_1*LAI))
    return ans
def PAR_GH():
    ans= tau_GH*n_GlobPar*I_Glob
    return ans
def PAR_FlrCan():
    ans=p_Flr*PAR_GH()*(1-p_Can)*math.exp(-K_1*LAI)*(1-math.exp(-K_2*LAI))
    return ans
def PAR_Can():
    ans = PAR_GHCan()+PAR_FlrCan()
    return ans
def CO2_Stom(CO2_Air):
    ans = n_CO2_AirStom*CO2_Air
    return ans
def Gam():
    ans = 1/LAI*c_Gam*T_Can+20*c_Gam*(1-1/LAI)
    return ans
def J():
    ans = (JPot()+alpha*PAR_Can()-((JPot()+alpha*PAR_Can())**2-4*Theta*JPot()*alpha*PAR_Can())**(1/2))/(2*Theta)    
    return ans
def PhoRate(CO2_Air):
    ans= J()*(CO2_Stom(CO2_Air)-Gam())/(4*(CO2_Stom(CO2_Air)+2*Gam()))
    return ans
def Res(CO2_Air):
    ans= PhoRate(CO2_Air)*(Gam()/CO2_Stom(CO2_Air))
    return ans 

def MC_AirCan(CO2_Air):
    ans = M_CH2O * h_CBuf() * (PhoRate(CO2_Air) - Res(CO2_Air))
    return ans

    # 19 


def h_CBuf():
    if C_Buf > C_MaxBuf:
        return 0
    else:
        return 1

    # 22


def P_calculate(Res, CO2_05, PMax,CO2Air):
    a = Res
    b = -(CO2Air + CO2_05 + Res * PMax)
    c = CO2Air * PMax

    if a == 0:
        return -c / b
    else:
        delta = b ** 2 - 4 * a * c
        res1 = (-b + math.sqrt(delta)) / (2 * a)
        res2 = (-b - math.sqrt(delta)) / (2 * a)

    # 24 doc them


def f_T(T):
    temp1 = (-Hd) / R * (1 / T_0 - S_entro / Hd)
    temp2 = (-Hd) / R * (1 / T - S_entro / Hd)

    ans = (1 + math.exp(temp1)) / (1 + math.exp(temp2))
    return ans

    # 25 doc them


def P_Max(T):
    return k_T_Arr(T) * f_T(T)

    # 27 doc them


def L_photon(K, m):
    ans = L0_27 * (1 - (K * math.exp(1) ** (-K * LAI)) / (1 - m))

    return ans

    # 28          J/mol||J/mol.K doc them


def k_T_Arr(T):
    ans = LAI * kT0_28 * math.exp(1) ** ((-Ha) / R * (T ** -1 - T_0 ** -1))

    return ans

    # 29 doc them


def P_Max_Mich(L, T):
    L_05 = P_Max(T) / 2
    ans = P_MLT * P_Max(T) * L
    ans /= (L + L_05)

    return ans

def AirDensity(h):
    return p_Air0*math.exp((g*h*M_Air)/(T_0*R))

# main
A_Flr = 1.4*(10**4) #
A_Roof = A_Flr*0.1  #
A_Side = A_Flr*0    #

# physics constant
g = 9.81  #
R=8.31447 #
M_Air = 28.96 #

LAI = 3  #

h_SideRoof = 0  #
h_Roof = 0.87  #
C_Buf = 1 #
h_Mean=4.2 #
h_Air=3.8 #

R = 8.314  #
T_0 = 293.15 #

S_entro = 710  # entropy
S_InsScr = 1  #

# concentration (mg/m3) = 0.0409 x concentration (ppm) x molecular weight
CO2_Air = 0.049*M_Air*692 # 
CO2_Top = CO2_Air   #
CO2_Out = 668  #  mg/m^3

o_Pad = 0
o_ExtCO2 = 720  #
o_VentForced=0

K_ThScr = 5 * 10**-4  #

T_Air = 19.9     +273.15   #
T_Top = T_Air+1 +273.15   #
T_Can = T_Air+1 +273.15   #
T_Out = 15.2    +273.15   #
T_MeanAir=(T_Air+T_Out)/2 +273.15 #
T_25 = T_0+25

p_Air0=1.2 #
p_MeanAir = AirDensity(h_Mean) #
p_Air = AirDensity(h_Air) # 
p_Top = p_Air+0.2 #
v_Wind = 4.6  #

c_leakage = 0.3 * (10**-4)  #

C_d = 0.35  #
C_w = 0.02  #

M_CH2O = 0.03       # #mg/mcmol
C_MaxBuf=20*(10**3) #
P_Blow=0.5*(10**6)  #
P=100

kT0_28 = 30
L0_27 = 3000
P_MLT = 7

capAir= 3.8 # h air
capTop= 0.4 # h green house - h air

U_ThScr = 0
U_Roof=1
U_Side=1
U_Pad = 0
U_Blow=0
U_ExtCO2=0
U_VentForced=0

n_Roof_Thr = 0.9 #vanthor
n_Side_Thr = 0.9 # bang voi n roof thr
n_HeatCO2 = 0.057  #
n_Roof  = 0.1 #
n_Side  =0 #

#canopy value in table 9.1 vanthor
H=22*(10**4)
S=710
K_1=0.7 
K_2=0.7
p_Can=0.07
tau_GH=0.78
n_GlobPar=2.3
I_Glob=1.8
p_Flr=0.5
Theta=0.7
c_Gam=1.7
n_CO2_AirStom=0.67
alpha = 0.385 
JMax = 210 
Ej=37*(10**3)
print("CO2_Air:",CO2_Air)
print("CO2_Top:",CO2_Top)

eulerval=[CO2_Air,CO2_Top]
rk4val=[CO2_Air,CO2_Top]
last = 50
min = 60
# h=300
# h=int(input("Step size:"))
# print(dx(CO2_Air,CO2_Top))
# print(euler(CO2_Air,CO2_Top,h))
# print(rk4(CO2_Air,CO2_Top,h))
##
# for j in range(5,last+5,5):
#     print("Phut thu:", j)
#     print("EULER:")
#     for i in range(0,(int)(min*60/h)):
#         eulerval= euler(eulerval[0],eulerval[1],h)
#     print(eulerval,sep='\n')
#     print("RK4:")
#     for i in range (0,(int)(min*60/h)):
#         rk4val=rk4(rk4val[0],rk4val[1],h)
#     print(rk4val,sep='\n')
##
#print("After 10 min")
#for i in range ((int)(5*60/h)):
#    printval= euler(printval[0],printval[1],h)
# #    print(printval,sep='\n')
# for i in range(0,300,1):
#     eulerval = euler(eulerval[0], eulerval[1],h)
#     rk4val=rk4(rk4val[0],rk4val[1],h)
    
# print(eulerval, rk4val)

def read_co2_rh():
    file = pd.read_csv("environment.csv")
    co2 = file.Co2[0]
    rh = file.Rh[0]
    return [co2,rh]

def read_co2_file():
    file = pd.read_csv("environment.csv")
    co2 = file.Co2
    return co2[0:576]

def plot(co2_rungekutta, co2_filecsv):
    # substraction = abs(co2_rungekutta-co2_filecsv)
    # plt.plot(substraction, label = 'co2_rungekutta-co2_filecsv')
    plt.plot(co2_rungekutta, label = 'Rungekutta')
    #plt.legend()
    plt.plot(co2_filecsv, label = 'FromFile')
    plt.legend()
    plt.axis([0,600,0,1000])
    plt.xlabel('Time every 5 mins')
    plt.ylabel('CO2')
    plt.title('Plot Co2 from file and by using rungekutta method')
    plt.show()

#main
co2_rk4 = [rk4val[0]]
for i in range(0,575,1):
    rk4val = rk4(rk4val[0],rk4val[1],300)
    co2_rk4 += [rk4val[0]]


co2_file = read_co2_file()
plot(co2_rk4,co2_file)