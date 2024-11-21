import sys
import os

project_dir = os.path.dirname(__file__)
if project_dir not in sys.path:
    sys.path.append(f'C:/Users/kylew/Documents/')

from SkewTProject import Bolton
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axisartist import Subplot
from matplotlib . ticker import FuncFormatter , Formatter
from mpl_toolkits . axisartist . grid_helper_curvelinear import GridHelperCurveLinear 
from SkewTProject.readsoundings import parse_SPC


C_to_K = 273.15
skew_slope = 40


def x_from_tp(T,p):
    x = T - skew_slope*np.log(p)
    return x

def y_from_p(p):
    y = - np.log(p)
    return y

def p_from_y(y):
    p = np.exp(-y)
    return p

def T_from_xp(x,p):
    T = x + skew_slope*np.log(p)
    return T

def to_thermo (x , y):
    """ Transform (x , y ) coordinates to T in degrees Celsius
    and p in mb . """
    p = p_from_y ( y )
    T_C = T_from_xp (x , p ) - C_to_K
    return T_C , p

def from_thermo ( T_C , p ):
    """ Transform T_C ( in degrees Celsius )
    and p ( in mb ) to (x , y ). """
    y = y_from_p ( p )
    x = x_from_tp (T_C + C_to_K , p)
    return x , y

# values along the bottom and left edges
p_bottom = 1050.0
p_top = 150
T_min = -40
T_max = 50

x_max, y_max = from_thermo(T_max, p_top)
x_min, y_min = from_thermo(T_min, p_bottom)

p_levels = np.arange(1000, 150-50, -50)
T_C_levels = np.arange(-80, 40+10, 10)
T_levels = T_C_levels + C_to_K
theta_levels = np.arange(-40, 100+10, 10) + C_to_K
theta_ep_levels = theta_levels.copy()
mixing_ratio =  np.asarray([.4, 1, 2, 3, 5, 8, 12, 16, 20])/1000

p_all = np.arange(p_bottom, p_top, -1)

y_p_levels = y_from_p(p_levels)
y_all_p = y_from_p(p_all)
x_T_levels = [x_from_tp(Ti, p_all) for Ti in T_levels]
x_thetas = [x_from_tp(Bolton.theta_dry(theta_i, p_all), p_all) for theta_i in theta_levels]
x_mixing_ratios = [x_from_tp(Bolton.mixing_ratio_line(p_all, mixing_ratio_i)+C_to_K, p_all) for mixing_ratio_i in mixing_ratio]
mesh_T , mesh_p = np.meshgrid(np.arange(-60.0,T_levels.max() - C_to_K +0.1, 0.1), p_all)
theta_ep_mesh = Bolton.theta_ep_field(mesh_T , mesh_p)

def Theta_E_Calc(T, p, p_0=1000):
    Rd = 287.04
    w_s = Bolton.sat_mixing_ratio(p, T)
    Tk = C_to_K + T
    Cl = 2336 + 1005.7
    Cdl = 1005.7 + w_s*Cl
    Rv = 461.5
    Lw = (2.5-(0.00237*T))*10**6
    rh = Bolton.RH(T,p,w_s)
    A = Tk*(p_0/p)**(Rd/Cdl)
    B = rh**((-w_s*Rv)/Cdl)
    C = np.exp((w_s*Lw)/(Cdl*Tk))
    return A*B*C

mesh_T , mesh_p = np.meshgrid(np.arange(-60.0,T_levels.max() - C_to_K +0.1, 0.1), p_all)
Theta_E_Mesh = Theta_E_Calc(mesh_T , mesh_p)

fig = plt . figure ( figsize =(12 ,12), dpi = 300)
skew_grid_helper = GridHelperCurveLinear((from_thermo, to_thermo))
ax = Subplot ( fig , 1 , 1 , 1, grid_helper = skew_grid_helper)
ax.set_xlabel('Temperature ( deg C )')
ax.set_ylabel('Pressure ( hPa )')
fig.add_subplot (ax)

for yi in y_p_levels:
    ax . plot (( x_min , x_max ) , ( yi , yi ) , color =(1.0 , 0.8 , 0.8))
for x_T in x_T_levels :
    ax . plot ( x_T , y_all_p , color =(1.0 , 0.5 , 0.5))
for x_theta in x_thetas :
    ax . plot ( x_theta , y_all_p , color =(1.0 , 0.7 , 0.7))
for x_mixing_ratio in x_mixing_ratios :
    good = p_all >= 600 # restrict mixing ratio lines to below 600 mb
    ax . plot ( x_mixing_ratio [ good ] , y_all_p [ good ] , color =(0.8 , .8 , 0.6))
n_moist = len ( theta_ep_levels )
moist_colors = ((0.6 ,0.9 ,0.7) ,)* n_moist
ax . contour ( x_from_tp ( mesh_T + C_to_K , mesh_p ) , y_from_p ( mesh_p ) ,
              theta_ep_mesh , theta_ep_levels , colors = moist_colors )
ax . axis (( x_min , x_max , y_min , y_max ))

n_moist = len ( theta_ep_levels )
moist_colors = ((0 ,0 ,0) ,)* n_moist
ax . contour ( x_from_tp ( mesh_T + C_to_K , mesh_p ) , y_from_p ( mesh_p ) ,
              Theta_E_Mesh , theta_levels , colors = moist_colors )
ax . axis (( x_min , x_max , y_min , y_max ))


def format_coord (x , y ):
    T , p = to_thermo (x , y )
    return " {0:5.1 f } C , {1:5.1 f } mb " . format ( float ( T ) , float ( p ))
ax . format_coord = format_coord

filepath = '/Users/kylew/Documents/SkewTProject/IAD.txt'

sounding_data = parse_SPC(filepath)

snd_T = sounding_data ['T']
# all temperature values , deg . C , should be in this range .


snd_P = sounding_data ['p']
good_P = ( snd_P > 150) & ( snd_P < 1010)
snd_P = snd_P[good_P]

good_T = ( snd_T > -100.0) & ( snd_T < 60.0)
snd_T = snd_T[good_P]


snd_Td = sounding_data ['Td']
good_Td = ( snd_Td > -100.0) & ( snd_Td < 60.0)
snd_Td = snd_Td[good_P]

x_snd_T = x_from_tp(snd_T+273.15, snd_P)
x_snd_Td = x_from_tp(snd_Td+273.15, snd_P)
y_snd_p = y_from_p(snd_P)

ax.plot( x_snd_Td , y_snd_p , linewidth =2 , color = 'g')
ax.plot( x_snd_T , y_snd_p , linewidth =2 , color = 'r')
plt.xlim(-50,50)
plt.show()
