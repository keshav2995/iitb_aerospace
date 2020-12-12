import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve

''' The given code takes m, t, p, alpha, M as inputs and returns 
    semi-wedge angles.
'''

class Teland:

    #Class variables
    m = 0.00     # value of maximum camber as a fraction of chord(0-1)'))
    t = 0.07     # value of maximum thickness as a fraction of chord(0-1)'))
    p = 0.5      # the location of maximum thickness as a fraction of chord(0-1)'))
    alpha = 2     # the value of angle of attack'))
    M1 = 3     # the value of free stream Mach Number'))
    gamma = 1.4   # value of adiabitc index of gas'))
    y_u = (t/2+m)
    y_l = -(t/2-m)
    theta_ul = np.arctan(y_u/p)*180/np.pi         # for upper left surface in degrees
    theta_ll = abs(np.arctan(y_l/p)*180/np.pi)     # for lower left surface in degrees
    theta_ur = np.arctan(y_u/(1-p))*180/np.pi      # for upper right surface in degrees
    theta_lr = abs(np.arctan(y_l/(1-p))*180/np.pi) # for lower right surface in degrees

    def __init__(self, x):
        #Instance Variables
        pass

    @classmethod
    def waveangle_u(cls,a):
        return np.tan((cls.theta_ul-cls.alpha)*np.pi/180) - 2*(np.cos(a)/np.sin(a))*(((cls.M1*np.sin(a))**2 - 1)/((cls.M1**2)*(1.4 + 2*(np.cos(a)**2) - 1)+2))

    @staticmethod
    def beta_values():
        B_u = fsolve(Teland.waveangle_u,0.00001)      #Value in radians  
        B1 = B_u*180/np.pi
        return B_u, B1  

    @classmethod
    def mu2_value(cls):
        Beta_u, Beta_1 = Teland.beta_values()
        Mn1 = cls.M1*np.sin(Beta_u)
        Mn2 = (((2/(cls.gamma-1))+Mn1**2)/(2*cls.gamma*Mn1**2/(cls.gamma-1)-1))**.5
        M2 = Mn2/np.sin((Beta_1-cls.theta_ul+cls.alpha)*np.pi/180)
        p2_p1 = 1+2*cls.gamma/(cls.gamma+1)*(Mn1**2-1)
        devtn = cls.theta_ul+cls.theta_ur
        mu1 = (((cls.gamma+1)/(cls.gamma-1))**.5*np.arctan(((cls.gamma-1)/(cls.gamma+1)*(M2**2-1))**.5)- np.arctan((M2**2-1)**.5))*180/np.pi
        mu2 = devtn+mu1
        return -mu2, M2, p2_p1

    @classmethod
    def pm_mach(cls, b):
        return Teland.mu2_value()[0]*np.pi/180 + (((cls.gamma+1)/(cls.gamma-1))**.5*(np.arctan(((cls.gamma-1)/(cls.gamma+1)*(b**2-1))**.5)) - np.arctan((b**2-1)**.5))


    @staticmethod
    def m3_value():
        M3 = fsolve(Teland.pm_mach, Teland.mu2_value()[1])
        return M3

    @classmethod
    def p3_p2_cal(cls):
        p3_p2 = ((1+.5*(cls.gamma-1)*Teland.mu2_value()[1]**2)/(1+.5*(cls.gamma-1)*Teland.m3_value()**2))**(cls.gamma/(cls.gamma-1))
        p3_p1 = p3_p2*Teland.mu2_value()[2]
        return p3_p1


    ''' Funtion computes the wave angle from above calculated wedge angles
        for lower surface.
    '''
    @classmethod
    def waveangle_l(cls, c):
        return np.tan((cls.theta_ll+cls.alpha)*np.pi/180) - 2*(np.cos(c)/np.sin(c))*(((cls.M1*np.sin(c))**2 - 1)/((cls.M1**2)*(1.4 + 2*(np.cos(c)**2) - 1)+2))
    
    
    @staticmethod
    def beta_values2():
        B_l = fsolve(Teland.waveangle_l,0.00001)      #Value in radians  
        B3 = B_l*180/np.pi                   # Value in degrees
        return B_l, B3

    @classmethod
    def interim_calc(cls):
        Mn11 = cls.M1*np.sin(Teland.beta_values2()[0])
        Mn4 = (((2/(cls.gamma-1))+Mn11**2)/(2*cls.gamma*Mn11**2/(cls.gamma-1)-1))**.5
        M4 = Mn4/np.sin((Teland.beta_values2()[1]-cls.theta_ll-cls.alpha)*np.pi/180)
        p4_p1 = 1+2*cls.gamma/(cls.gamma+1)*(Mn11**2-1)
        devtn2 = cls.theta_ll+cls.theta_lr
        mu4 = (((cls.gamma+1)/(cls.gamma-1))**.5*np.arctan(((cls.gamma-1)/(cls.gamma+1)*(M4**2-1))**.5)- np.arctan((M4**2-1)**.5))*180/np.pi
        mu5 = devtn2+mu4
        return -mu5, M4, p4_p1

    @classmethod
    def pm_mach1(cls,d):
        return Teland.interim_calc()[0]*np.pi/180 + (((cls.gamma+1)/(cls.gamma-1))**.5*(np.arctan(((cls.gamma-1)/(cls.gamma+1)*(d**2-1))**.5)) - np.arctan((d**2-1)**.5))

    @staticmethod
    def m5_value():
        M5 = fsolve(Teland.pm_mach1, Teland.interim_calc()[1])
        return M5

    @classmethod
    def final_calc(cls):
        p5_p4 = ((1+.5*(cls.gamma-1)*Teland.interim_calc()[1]**2)/(1+.5*(cls.gamma-1)*Teland.m5_value()**2))**(cls.gamma/(cls.gamma-1))
        p5_p1 = p5_p4*Teland.interim_calc()[2]

        Cn = ((Teland.interim_calc()[2]*(np.cos(cls.theta_ll*np.pi/180))-Teland.mu2_value()[2]*(np.cos(cls.theta_ul*np.pi/180)))*(cls.p)+
        (p5_p1*(np.cos(cls.theta_lr*np.pi/180))-Teland.p3_p2_cal()*(np.cos(cls.theta_ur*np.pi/180)))*(1-cls.p))/(.5*cls.M1**2*cls.gamma)
        
        Ca = ((Teland.interim_calc()[2]*(np.sin(cls.theta_ll*np.pi/180))+Teland.mu2_value()[2]*(np.sin(cls.theta_ul*np.pi/180)))*(cls.p)-
        (p5_p1*(np.sin(cls.theta_lr*np.pi/180))+Teland.p3_p2_cal()*(np.sin(cls.theta_ur*np.pi/180)))*(1-cls.p))/(.5*cls.M1**2*cls.gamma)
        
        alpha1 = cls.alpha * np.pi/180

        Cl = Cn*np.cos(alpha1) - Ca*np.sin(alpha1)
        Cd = Cn*np.sin(alpha1) + Ca*np.cos(alpha1)
        return Cl, Cd

    @classmethod
    def ackeret_solution(cls):
        slopeuf = np.tan((cls.theta_ul-cls.alpha)*np.pi/180)
        slopeua = np.tan((-cls.alpha-cls.theta_ur)*np.pi/180)
        slopelf = np.tan((-cls.alpha-cls.theta_ll)*np.pi/180)
        slopela = np.tan((-cls.alpha+cls.theta_lr)*np.pi/180)
     
        cp2 = 2*slopeuf/np.sqrt(cls.M1**2 -1)
        cp3 = 2*slopeua/np.sqrt(cls.M1**2 -1)
        cp4 = -2*slopelf/np.sqrt(cls.M1**2 -1)
        cp5 = -2*slopela/np.sqrt(cls.M1**2 -1)

        cn = (cp4*(np.cos(cls.theta_ll*np.pi/180)) - cp2*(np.cos(cls.theta_ul*np.pi/180)))*cls.p + (cp5*(np.cos(cls.theta_lr*np.pi/180)) - cp3*(np.cos(cls.theta_ur*np.pi/180)))*(1-cls.p)
        ca = (cp4*(np.sin(cls.theta_ll*np.pi/180)) + cp2*(np.sin(cls.theta_ul*np.pi/180)))*cls.p - (cp5*(np.sin(cls.theta_lr*np.pi/180)) + cp3*(np.sin(cls.theta_ur*np.pi/180)))*(1-cls.p)

        alpha1= cls.alpha*np.pi/180
            
        Cl = cn*np.cos(alpha1)-ca*np.sin(alpha1)
        Cd = cn*np.sin(alpha1)+ca*np.cos(alpha1)
        return Cl, Cd

def main():

    C_lift, C_drag = Teland.final_calc()
    print("(i)Sectional lift coefficient: {} \n (ii) Sectional wave drag coefficient: {}".format(C_lift[0], C_drag[0]))

    C_lift, C_drag = Teland.ackeret_solution()
    print("(i)Sectional lift coefficient by AT: {} \n (ii) Sectional wave drag coefficient by AT: {}".format(C_lift, C_drag))

if __name__ == "__main__":
    main()
