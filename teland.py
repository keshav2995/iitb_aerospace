import numpy as np
import pandas as pd
import plotly.express as px

class Teland:
    
    naca = '1312'
    
    #p is loc_max_camber = 3/10
    #m is max_camber = 1/100
    #x is in range of [0,1]  

    def __init__(self):
        self.loc_max_camber = int(self.naca[1])/10
        self.max_camber = int(self.naca[0])/100
        self.x = 0
        self.plots = []

    def z_c(self):
        for self.x in np.linspace(0.0,1.0,100):
            if self.x <= self.loc_max_camber: 
                self.Z_c = (self.max_camber)*(2*self.loc_max_camber*self.x - self.x**2)/(self.loc_max_camber)**2
                self.plots.append([self.Z_c,self.x])
            else:
                self.Z_c = (self.max_camber)*(1-2*self.loc_max_camber+2*self.loc_max_camber*self.x - self.x**2)/(1-self.loc_max_camber)**2
                self.plots.append([self.Z_c,self.x]) 
        return self.plots


def main():

    tl = Teland()
    plots = tl.z_c()
    df = pd.DataFrame(plots, columns=['Camber Function(Zc) w.r.t x', 'x'])
    fig = px.line(df, x="x", y="Camber Function(Zc) w.r.t x", title='Teland ka Graph')
    fig.show()


if __name__ == "__main__":
    main()   

