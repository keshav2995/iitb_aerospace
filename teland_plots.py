import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve
import pandas as pd 
import seaborn as sns 
import matplotlib.pyplot as plt 
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

''' The given code takes m, t, p, alpha, M as inputs and returns 
    semi-wedge angles.
'''

class Teland:

    gamma = 1.4   # value of adiabitc index of gas'))

    def __init__(self, m, t, p, alpha, M1, var_plots=False):
        if not var_plots:
            self.m = m     # value of maximum camber as a fraction of chord(0-1)'))
            self.t = t     # value of maximum thickness as a fraction of chord(0-1)'))
            self.p = p      # the location of maximum thickness as a fraction of chord(0-1)'))
            self.alpha = alpha     # the value of angle of attack'))
            self.M1 = M1     # the value of free stream Mach Number'))

        else:
            self.M1 = np.linspace(int(M1[0]),int(M1[1]),10)
            self.t = np.linspace(float(t[0]),float(t[1]),10)
            self.p = np.linspace(float(p[0]),float(p[1]),10)
            self.alpha = np.linspace(int(alpha[0]),int(alpha[1]),10)
            self.m = np.linspace(float(m[0]),float(m[1]),10)

    @classmethod
    def from_user_input(self):
        return self(
            input('Enter range of M separated by comma: ').split(','), #1,4
            input('Enter range of t separated by comma: ').split(','), #0.03,0.06
            input('Enter range of p separated by comma:: ').split(','), #0.1,0.5
            input('Enter range of alpha separated by comma:: ').split(','), #0,5
            input('Enter range of m separated by comma:: ').split(','), #-0.01,0.01
            input('Enter "True" for plots:')
        )

    def initial_variables(self):
        self.y_u = (self.t/2+self.m)
        self.y_l = -(self.t/2-self.m)
        self.theta_ul = np.arctan(self.y_u/self.p)*180/np.pi         # for upper left surface in degrees
        self.theta_ll = abs(np.arctan(self.y_l/self.p)*180/np.pi)     # for lower left surface in degrees
        self.theta_ur = np.arctan(self.y_u/(1-self.p))*180/np.pi      # for upper right surface in degrees
        self.theta_lr = abs(np.arctan(self.y_l/(1-self.p))*180/np.pi) # for lower right surface in degrees
        

    def waveangle_u(self,a):
        return np.tan((self.theta_ul-self.alpha)*np.pi/180) - 2*(np.cos(a)/np.sin(a))*(((self.M1*np.sin(a))**2 - 1)/((self.M1**2)*(1.4 + 2*(np.cos(a)**2) - 1)+2))

    
    def beta_values(self):
        self.B_u = fsolve(self.waveangle_u,0.00001)      #Value in radians  
        self.B1 = self.B_u*180/np.pi
        return self.B_u, self.B1  

    def mu2_value(self):
        self.Beta_u, self.Beta_1 = self.beta_values()
        self.Mn1 = self.M1*np.sin(self.Beta_u)
        self.Mn2 = (((2/(self.gamma-1))+self.Mn1**2)/(2*self.gamma*self.Mn1**2/(self.gamma-1)-1))**.5
        self.M2 = self.Mn2/np.sin((self.Beta_1-self.theta_ul+self.alpha)*np.pi/180)
        self.p2_p1 = 1+2*self.gamma/(self.gamma+1)*(self.Mn1**2-1)
        self.devtn = self.theta_ul+self.theta_ur
        self.mu1 = (((self.gamma+1)/(self.gamma-1))**.5*np.arctan(((self.gamma-1)/(self.gamma+1)*(self.M2**2-1))**.5)- np.arctan((self.M2**2-1)**.5))*180/np.pi
        self.mu2 = self.devtn+self.mu1

    def pm_mach(self, b):
        return -self.mu2*np.pi/180 + (((self.gamma+1)/(self.gamma-1))**.5*(np.arctan(((self.gamma-1)/(self.gamma+1)*(b**2-1))**.5)) - np.arctan((b**2-1)**.5))

    def m3_value(self):
        self.M3 = fsolve(self.pm_mach1, self.M2)

    def p3_p2_cal(self):
        self.p3_p2 = ((1+.5*(self.gamma-1)*self.M2**2)/(1+.5*(self.gamma-1)*self.M3**2))**(self.gamma/(self.gamma-1))
        self.p3_p1 = self.p3_p2*self.p2_p1

    ''' Funtion computes the wave angle from above calculated wedge angles
        for lower surface.
    '''
    def waveangle_l(self, c):
        return np.tan((self.theta_ll+self.alpha)*np.pi/180) - 2*(np.cos(c)/np.sin(c))*(((self.M1*np.sin(c))**2 - 1)/((self.M1**2)*(1.4 + 2*(np.cos(c)**2) - 1)+2))
    
    def beta_values2(self):
        self.B_l = fsolve(self.waveangle_l,0.00001)      #Value in radians  
        self.B3 = self.B_l*180/np.pi                   # Value in degrees

    def interim_calc(self):
        self.Mn11 = self.M1*np.sin(self.B_l)
        self.Mn4 = (((2/(self.gamma-1))+self.Mn11**2)/(2*self.gamma*self.Mn11**2/(self.gamma-1)-1))**.5
        self.M4 = self.Mn4/np.sin((self.B3-self.theta_ll-self.alpha)*np.pi/180)
        self.p4_p1 = 1+2*self.gamma/(self.gamma+1)*(self.Mn11**2-1)
        self.devtn2 = self.theta_ll+self.theta_lr
        self.mu4 = (((self.gamma+1)/(self.gamma-1))**.5*np.arctan(((self.gamma-1)/(self.gamma+1)*(self.M4**2-1))**.5)- np.arctan((self.M4**2-1)**.5))*180/np.pi
        self.mu5 = self.devtn2+self.mu4

    def pm_mach1(self,d):
        return -self.mu5*np.pi/180 + (((self.gamma+1)/(self.gamma-1))**.5*(np.arctan(((self.gamma-1)/(self.gamma+1)*(d**2-1))**.5)) - np.arctan((d**2-1)**.5))

    def m5_value(self):
        self.M5 = fsolve(self.pm_mach1, self.M4)

    def final_calc(self):
        self.p5_p4 = ((1+.5*(self.gamma-1)*self.M4**2)/(1+.5*(self.gamma-1)*self.M5**2))**(self.gamma/(self.gamma-1))
        self.p5_p1 = self.p5_p4*self.p4_p1

        self.Cn = ((self.p4_p1*(np.cos(self.theta_ll*np.pi/180))-self.p2_p1*(np.cos(self.theta_ul*np.pi/180)))*(self.p)+
        (self.p5_p1*(np.cos(self.theta_lr*np.pi/180))-self.p3_p1*(np.cos(self.theta_ur*np.pi/180)))*(1-self.p))/(.5*self.M1**2*self.gamma)
        
        self.Ca = ((self.p4_p1*(np.sin(self.theta_ll*np.pi/180))+self.p2_p1*(np.sin(self.theta_ul*np.pi/180)))*(self.p)-
        (self.p5_p1*(np.sin(self.theta_lr*np.pi/180))+self.p3_p1*(np.sin(self.theta_ur*np.pi/180)))*(1-self.p))/(.5*self.M1**2*self.gamma)
        
        self.alpha1 = self.alpha * np.pi/180

        self.Cl = self.Cn*np.cos(self.alpha1) - self.Ca*np.sin(self.alpha1)
        self.Cd = self.Cn*np.sin(self.alpha1) + self.Ca*np.cos(self.alpha1)

    # Solution using Ackeret Theory
    def ackeret_solution(self):
        self.slopeuf = np.tan((self.theta_ul-self.alpha)*np.pi/180)
        self.slopeua = np.tan((-self.alpha-self.theta_ur)*np.pi/180)
        self.slopelf = np.tan((-self.alpha-self.theta_ll)*np.pi/180)
        self.slopela = np.tan((-self.alpha+self.theta_lr)*np.pi/180)
     
        self.cp2 = 2*self.slopeuf/np.sqrt(self.M1**2 -1)
        self.cp3 = 2*self.slopeua/np.sqrt(self.M1**2 -1)
        self.cp4 = -2*self.slopelf/np.sqrt(self.M1**2 -1)
        self.cp5 = -2*self.slopela/np.sqrt(self.M1**2 -1)

        self.cn = (self.cp4*(np.cos(self.theta_ll*np.pi/180)) - self.cp2*(np.cos(self.theta_ul*np.pi/180)))*self.p + (self.cp5*(np.cos(self.theta_lr*np.pi/180)) - self.cp3*(np.cos(self.theta_ur*np.pi/180)))*(1-self.p)
        self.ca = (self.cp4*(np.sin(self.theta_ll*np.pi/180)) + self.cp2*(np.sin(self.theta_ul*np.pi/180)))*self.p - (self.cp5*(np.sin(self.theta_lr*np.pi/180)) + self.cp3*(np.sin(self.theta_ur*np.pi/180)))*(1-self.p)

        self.alpha1= self.alpha*np.pi/180
            
        self.Cl = self.cn*np.cos(self.alpha1)-self.ca*np.sin(self.alpha1)
        self.Cd = self.cn*np.sin(self.alpha1)+self.ca*np.cos(self.alpha1)


    def visualization(self):
        #-------*******todo*******--------#
        for i in self.var_M:
            cl_,cd_ = Teland.final_calc()
        fig = make_subplots(rows=2, cols=2,
                            specs=[[{"secondary_y": True}, {"secondary_y": True}],
                                [{"secondary_y": True}, {"secondary_y": True}]])

        # 1st left
        fig.add_trace(
            go.Scatter(x=self.var_m, y=Teland.final_calc(), name="CL/CD vs m for var M"),
            row=1, col=1, secondary_y=False,
        )

        fig.add_trace(
            go.Scatter(x=df['Time'], y=df['BB'], name="BB"),
            row=1, col=1, secondary_y=True,
        )

        # 1rd right
        fig.add_trace(
            go.Scatter(x=df['Time'], y=df['Price'], name="Price"),
            row=1, col=2, secondary_y=False,
        )

        fig.add_trace(
            go.Scatter(x=df['Time'], y=df['MACD'], name="MACD"),
            row=1, col=2, secondary_y=True,
        )

        # Bottom left
        fig.add_trace(
            go.Scatter(x=df['Time'], y=df['Price'], name="Price"),
            row=2, col=1, secondary_y=False,
        )

        fig.add_trace(
            go.Scatter(x=df['Time'], y=df['EMA'], name="EMA"),
            row=2, col=1, secondary_y=True,
        )

        # Bottom right
        fig.add_trace(
            go.Scatter(x=df['Time'], y=df['Price'], name="Price"),
            row=2, col=2, secondary_y=False,
        )

        fig.add_trace(
            go.Scatter(x=df['Time'], y=df['MA'], name="MA"),
            row=2, col=2, secondary_y=True,
        )


        fig.show()



def main():
    t1 = Teland(0.00,0.07,0.5,2,3)
    C_lift, C_drag = Teland.final_calc()
    print("(i)Sectional lift coefficient: {} \n (ii) Sectional wave drag coefficient: {}".format(C_lift[0], C_drag[0]))

    C_lift_at, C_drag_at = Teland.ackeret_solution()
    print("(i)Sectional lift coefficient by AT: {} \n (ii) Sectional wave drag coefficient by AT: {}".format(C_lift_at, C_drag_at))

    #Visualization
    # cl_ci = []
    # for i in np.linspace(0.03,0.06,10):
    #     C_lift, C_drag = Teland.final_calc()
    #     ratio = C_lift/C_drag
    #     cl_ci.append(ratio)


if __name__ == "__main__":
    main()
