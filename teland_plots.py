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

    def __init__(self, m, t, p, alpha, M1):
        self.m = m     # value of maximum camber as a fraction of chord(0-1)'))
        self.t = t     # value of maximum thickness as a fraction of chord(0-1)'))
        self.p = p      # the location of maximum thickness as a fraction of chord(0-1)'))
        self.alpha = alpha     # the value of angle of attack'))
        self.M1 = M1     # the value of free stream Mach Number'))

        # else:
        #     self.m = np.linspace(-0.01,0.01,25)
        #     self.t = np.linspace(0.03,0.06,5)
        #     self.p = np.linspace(0.1,0.5,5)
        #     self.alpha = np.linspace(0,5,5)
        #     self.M1 = np.linspace(1,4,5)
            # self.m = np.linspace(float(m[0]),float(m[1]),10)
            # self.t = np.linspace(float(t[0]),float(t[1]),10)
            # self.p = np.linspace(float(p[0]),float(p[1]),10)
            # self.alpha = np.linspace(int(alpha[0]),int(alpha[1]),10)
            # self.M1 = np.linspace(int(M1[0]),int(M1[1]),10)

    # @classmethod
    # def from_user_input(self):
    #     return self(
    #         input('Enter range of m separated by comma:: ').split(','), #-0.01,0.01
    #         input('Enter range of t separated by comma: ').split(','), #0.03,0.06
    #         input('Enter range of p separated by comma:: ').split(','), #0.1,0.5
    #         input('Enter range of alpha separated by comma:: ').split(','), #0,5
    #         input('Enter range of M separated by comma: ').split(','), #1,4
    #         input('Enter "True" for plots:')
    #     )

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
        self.M3 = fsolve(self.pm_mach, self.M2)

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
        return self.B_l, self.B3

    def interim_calc(self):
        self.Beta_l, self.Beta_3 = self.beta_values2()
        self.Mn11 = self.M1*np.sin(self.Beta_l)
        self.Mn4 = (((2/(self.gamma-1))+self.Mn11**2)/(2*self.gamma*self.Mn11**2/(self.gamma-1)-1))**.5
        self.M4 = self.Mn4/np.sin((self.Beta_3-self.theta_ll-self.alpha)*np.pi/180)
        self.p4_p1 = 1+2*self.gamma/(self.gamma+1)*(self.Mn11**2-1)
        self.devtn2 = self.theta_ll+self.theta_lr
        self.mu4 = (((self.gamma+1)/(self.gamma-1))**.5*np.arctan(((self.gamma-1)/(self.gamma+1)*(self.M4**2-1))**.5)- np.arctan((self.M4**2-1)**.5))*180/np.pi
        self.mu5 = self.devtn2+self.mu4

    def pm_mach1(self,d):
        return -self.mu5*np.pi/180 + (((self.gamma+1)/(self.gamma-1))**.5*(np.arctan(((self.gamma-1)/(self.gamma+1)*(d**2-1))**.5)) - np.arctan((d**2-1)**.5))

    def m5_value(self):
        self.M5 = fsolve(self.pm_mach1, self.M4)

    def runall(self):
        self.initial_variables()
        self.beta_values()
        self.mu2_value()
        self.m3_value()
        self.p3_p2_cal()
        self.beta_values2()
        self.interim_calc()
        self.m5_value()

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
        return self.Cl, self.Cd


    # Solution using Ackeret Theory
    def ackeret_solution(self):
        self.initial_variables()
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
        return self.Cl, self.Cd

    @staticmethod
    def visualization(df_t, df_p, df_alpha, df_M):
        # Create traces

        fig = make_subplots(rows=4, cols=2,
                    subplot_titles=("Cl/Cd vs m (for M) [SE theory]", 
                    "Cl/Cd vs m (for M) [Ackeret theory]", 
                    "Cl/Cd vs m (for t) [SE theory]", 
                    "Cl/Cd vs m (for t) [Ackeret theory]",
                    "Cl/Cd vs m (for p) [SE theory]",
                    "Cl/Cd vs m (for p) [Ackeret theory]",
                    "Cl/Cd vs m (for alpha) [SE theory]",
                    "Cl/Cd vs m (for alpha) [Ackeret theory]"
                    ),
                    specs=[[{"secondary_y": True}, {"secondary_y": True}],
                           [{"secondary_y": True}, {"secondary_y": True}],
                           [{"secondary_y": True}, {"secondary_y": True}],
                           [{"secondary_y": True}, {"secondary_y": True}]])

        fig.add_trace(go.Scatter(x=df_M.iloc[:,0], y=df_M.iloc[:,1],
                            mode='lines',
                            name=df_M.columns[1]),
                            row=1, col=1, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_M.iloc[:,0], y=df_M.iloc[:,3],
                            mode='lines',
                            name=df_M.columns[3]),
                            row=1, col=1, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_M.iloc[:,0], y=df_M.iloc[:,5],
                            mode='lines', 
                            name=df_M.columns[5]),
                            row=1, col=1, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_M.iloc[:,0], y=df_M.iloc[:,7],
                            mode='lines', 
                            name=df_M.columns[7]),
                            row=1, col=1, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_M.iloc[:,0], y=df_M.iloc[:,9],
                            mode='lines', 
                            name=df_M.columns[9]),
                            row=1, col=1, secondary_y=False,)                         

        fig.add_trace(go.Scatter(x=df_M.iloc[:,0], y=df_M.iloc[:,2],
                            mode='lines',
                            name=df_M.columns[2]),
                            row=1, col=2, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_M.iloc[:,0], y=df_M.iloc[:,4],
                            mode='lines',
                            name=df_M.columns[4]),
                            row=1, col=2, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_M.iloc[:,0], y=df_M.iloc[:,6],
                            mode='lines', 
                            name=df_M.columns[6]),
                            row=1, col=2, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_M.iloc[:,0], y=df_M.iloc[:,8],
                            mode='lines', 
                            name=df_M.columns[8]),
                            row=1, col=2, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_M.iloc[:,0], y=df_M.iloc[:,10],
                            mode='lines', 
                            name=df_M.columns[10]),
                            row=1, col=2, secondary_y=False,)                        

        fig.add_trace(go.Scatter(x=df_t.iloc[:,0], y=df_t.iloc[:,1],
                            mode='lines',
                            name=df_t.columns[1]),
                            row=2, col=1, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_t.iloc[:,0], y=df_t.iloc[:,3],
                            mode='lines',
                            name=df_t.columns[3]),
                            row=2, col=1, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_t.iloc[:,0], y=df_t.iloc[:,5],
                            mode='lines', 
                            name=df_t.columns[5]),
                            row=2, col=1, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_t.iloc[:,0], y=df_t.iloc[:,7],
                            mode='lines', 
                            name=df_t.columns[7]),
                            row=2, col=1, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_t.iloc[:,0], y=df_t.iloc[:,9],
                            mode='lines', 
                            name=df_t.columns[9]),
                            row=2, col=1, secondary_y=False,)                        

        fig.add_trace(go.Scatter(x=df_t.iloc[:,0], y=df_t.iloc[:,2],
                            mode='lines',
                            name=df_t.columns[2]),
                            row=2, col=2, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_t.iloc[:,0], y=df_t.iloc[:,4],
                            mode='lines',
                            name=df_t.columns[4]),
                            row=2, col=2, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_t.iloc[:,0], y=df_t.iloc[:,6],
                            mode='lines', 
                            name=df_t.columns[6]),
                            row=2, col=2, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_t.iloc[:,0], y=df_t.iloc[:,8],
                            mode='lines', 
                            name=df_t.columns[8]),
                            row=2, col=2, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_t.iloc[:,0], y=df_t.iloc[:,10],
                            mode='lines', 
                            name=df_t.columns[10]),
                            row=2, col=2, secondary_y=False,)                        

        fig.add_trace(go.Scatter(x=df_p.iloc[:,0], y=df_p.iloc[:,1],
                            mode='lines',
                            name=df_p.columns[1]),
                            row=3, col=1, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_p.iloc[:,0], y=df_p.iloc[:,3],
                            mode='lines',
                            name=df_p.columns[3]),
                            row=3, col=1, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_p.iloc[:,0], y=df_p.iloc[:,5],
                            mode='lines', 
                            name=df_p.columns[5]),
                            row=3, col=1, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_p.iloc[:,0], y=df_p.iloc[:,7],
                            mode='lines', 
                            name=df_p.columns[7]),
                            row=3, col=1, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_p.iloc[:,0], y=df_p.iloc[:,9],
                            mode='lines', 
                            name=df_p.columns[9]),
                            row=3, col=1, secondary_y=False,)                        

        fig.add_trace(go.Scatter(x=df_p.iloc[:,0], y=df_p.iloc[:,2],
                            mode='lines',
                            name=df_p.columns[2]),
                            row=3, col=2, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_p.iloc[:,0], y=df_p.iloc[:,4],
                            mode='lines',
                            name=df_p.columns[4]),
                            row=3, col=2, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_p.iloc[:,0], y=df_p.iloc[:,6],
                            mode='lines', 
                            name=df_p.columns[6]),
                            row=3, col=2, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_p.iloc[:,0], y=df_p.iloc[:,8],
                            mode='lines', 
                            name=df_p.columns[8]),
                            row=3, col=2, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_p.iloc[:,0], y=df_p.iloc[:,10],
                            mode='lines', 
                            name=df_p.columns[10]),
                            row=3, col=2, secondary_y=False,)                        

        fig.add_trace(go.Scatter(x=df_alpha.iloc[:,0], y=df_alpha.iloc[:,1],
                            mode='lines',
                            name=df_alpha.columns[1]),
                            row=4, col=1, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_alpha.iloc[:,0], y=df_alpha.iloc[:,3],
                            mode='lines',
                            name=df_alpha.columns[3]),
                            row=4, col=1, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_alpha.iloc[:,0], y=df_alpha.iloc[:,5],
                            mode='lines', 
                            name=df_alpha.columns[5]),
                            row=4, col=1, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_alpha.iloc[:,0], y=df_alpha.iloc[:,7],
                            mode='lines', 
                            name=df_alpha.columns[7]),
                            row=4, col=1, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_alpha.iloc[:,0], y=df_alpha.iloc[:,9],
                            mode='lines', 
                            name=df_alpha.columns[9]),
                            row=4, col=1, secondary_y=False,)                         

        fig.add_trace(go.Scatter(x=df_alpha.iloc[:,0], y=df_alpha.iloc[:,2],
                            mode='lines',
                            name=df_alpha.columns[2]),
                            row=4, col=2, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_alpha.iloc[:,0], y=df_alpha.iloc[:,4],
                            mode='lines',
                            name=df_alpha.columns[4]),
                            row=4, col=2, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_alpha.iloc[:,0], y=df_alpha.iloc[:,6],
                            mode='lines', 
                            name=df_alpha.columns[6]),
                            row=4, col=2, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_alpha.iloc[:,0], y=df_alpha.iloc[:,8],
                            mode='lines', 
                            name=df_alpha.columns[8]),
                            row=4, col=2, secondary_y=False,)

        fig.add_trace(go.Scatter(x=df_alpha.iloc[:,0], y=df_alpha.iloc[:,10],
                            mode='lines', 
                            name=df_alpha.columns[10]),
                            row=4, col=2, secondary_y=False,)                        

        # Update xaxis properties
        fig.update_xaxes(title_text="Values of m", row=1, col=1)
        fig.update_xaxes(title_text="Values of m", row=1, col=2)
        fig.update_xaxes(title_text="Values of m", row=2, col=1)
        fig.update_xaxes(title_text="Values of m", row=2, col=2)
        fig.update_xaxes(title_text="Values of m", row=3, col=1)
        fig.update_xaxes(title_text="Values of m", row=3, col=2)
        fig.update_xaxes(title_text="Values of m", row=4, col=1)
        fig.update_xaxes(title_text="Values of m", row=4, col=2)

        # Update yaxis properties
        fig.update_yaxes(title_text="Cl/Cd ratio", row=1, col=1)
        fig.update_yaxes(title_text="Cl/Cd ratio", row=1, col=2)
        fig.update_yaxes(title_text="Cl/Cd ratio", row=2, col=1)
        fig.update_yaxes(title_text="Cl/Cd ratio", row=2, col=2)
        fig.update_yaxes(title_text="Cl/Cd ratio", row=3, col=1)
        fig.update_yaxes(title_text="Cl/Cd ratio", row=3, col=2)
        fig.update_yaxes(title_text="Cl/Cd ratio", row=4, col=1)
        fig.update_yaxes(title_text="Cl/Cd ratio", row=4, col=2)

        fig.layout.update(height=1080, width=1280, showlegend=False, hovermode='closest')

        fig.show()



def main():
    t1 = Teland(0.00,0.07,0.5,2,3)
    t1.runall()
    C_lift, C_drag = t1.final_calc()
    print("(i)Sectional lift coefficient: {} \n (ii) Sectional wave drag coefficient: {}".format(C_lift[0], C_drag[0]))

    t2 = Teland(0.00,0.07,0.5,2,3)
    C_lift_at, C_drag_at = t2.ackeret_solution()
    print("(i)Sectional lift coefficient by AT: {} \n (ii) Sectional wave drag coefficient by AT: {}".format(C_lift_at, C_drag_at))

    seq = [1,2,3,4]
    names = ['t','p','alpha','M']
    #for cl/cd vs m (given val)
    def get_plot_data(*args):  #m, t, p, alpha, M1
        df_x = pd.DataFrame()
        arg = args
        df_x['m value'] = arg[0]
        global i_v
        i_ = i_v
        for x in arg[i_]:
            ratio_cl_cd_set = []
            ratio_cl_cd_at = []
            for y in arg[0]: #---> m
                if i_ ==1: 
                    t3 = Teland(m = y, t = x, p = arg[2], alpha = arg[3], M1 = arg[4]) 
                    t4 = Teland(y,x,arg[2],arg[3],arg[4])
                elif i_==2:
                    t3 = Teland(m = y, t = arg[1], p = x, alpha = arg[3], M1 = arg[4])
                    t4 = Teland(y,arg[1],x,arg[3],arg[4])
                elif i_==3:
                    t3 = Teland(m = y, t = arg[1], p = arg[2], alpha = x, M1 = arg[4])
                    t4 = Teland(y,arg[1],arg[2],x,arg[4])
                else:
                    t3 = Teland(m = y, t = arg[1], p = arg[2], alpha = arg[3], M1 = x)
                    t4 = Teland(y,arg[1],arg[2],arg[3],x)
                t3.runall()
                C_l, C_d = t3.final_calc()
                ratio = C_l[0]/C_d[0] 
                ratio_cl_cd_set.append(ratio)

                C_l_, C_d_ = t4.ackeret_solution()
                ratio_ = C_l_/C_d_
                ratio_cl_cd_at.append(ratio_)

            df_x['Cl/Cd ratio for {}={} (set)'.format(names[i_-1],x)] = ratio_cl_cd_set
            df_x['Cl/Cd ratio for {}={} (at)'.format(names[i_-1],x)] = ratio_cl_cd_at
            
        i_v += 1

        return df_x

    for i in seq:
        if i==4:
            df_M = get_plot_data(np.linspace(-0.01,0.01,25), 0.07, 0.5, 2, np.linspace(2,5,5))
        elif i==1:
            df_t = get_plot_data(np.linspace(-0.01,0.01,25), np.linspace(0.03,0.06,5), 0.5, 2, 3)
        elif i==2:
            df_p = get_plot_data(np.linspace(-0.01,0.01,25), 0.07, np.linspace(0.1,0.5,5), 2, 3)
        elif i==3:
            df_alpha = get_plot_data(np.linspace(-0.01,0.01,25), 0.07, 0.5, np.linspace(0,5,5), 3)


    Teland.visualization(df_t, df_p, df_alpha, df_M)

if __name__ == "__main__":
    i_v = 1
    main()
