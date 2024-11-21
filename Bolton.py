import numpy as np

C_to_K = 273.15
c_p_dry = 1005.7
c_V_dry = 718.66
eps = 0.622
k_dry = 0.2854

def sat_vapor_pressure(T):
    if np.any(T) > 200:
        print('Uh oh')
    e_s = 6.112 * np.exp((17.67*T)/(T+243.5))
    return e_s

def sat_vapor_temperature(e_s):
    T = (243.5*np.log(e_s)-440.8)/(19.48-np.log(e_s))
    return T

def sat_mixing_ratio(p, T):
    #p in mb and T in C
    e_s = sat_vapor_pressure(T)
    w = eps*((e_s)/(p-e_s))
    return w

def mixing_ratio_line(p, w_s):
    #pressure in mb and saturation mixing ratio in kg/kg
    e_s = (p*w_s)/(eps+w_s)
    T = sat_vapor_temperature(e_s)
    return T

def RH(T, p, w):
    #Temperature in C, pressure in mb, mixing ratio in kg/kg
    e_s = sat_vapor_pressure(T)
    e = (p*w)/(eps+w)
    return e/e_s

def T_LCL(T, RH):
    #T in kelvin and RH in %
    T_LCL = (1)/(((1)/(T-55))-(np.log(RH/100))/2840)
    return T_LCL

def theta_dry(theta, p, p_0=1000):
    #Theta in K, Pressure in mb
    e = 0
    T = theta*(((p-e)/p_0)**k_dry)
    return T

def pseudoeq_potential_T(T,p,w,p_0=1000):
    #Temp in C, Pressure in mb, mixing ratio in kg/kg
    Tk = T + C_to_K
    rh = RH(T,p,w)
    Tl = T_LCL(Tk, rh)
    A = Tk*(p_0/p)**(0.2854*(1-(0.28*10**-3)*w))
    B = np.exp(((3.376/Tl)-0.00254) * w*(1+(0.81*10**-3)*w))
    return A*B

def theta_ep_field(T, p, p_0=1000.0):
    w_s = sat_mixing_ratio(p, T)*1000
    theta_ep = pseudoeq_potential_T(T, p, w_s, p_0)
    return theta_ep


