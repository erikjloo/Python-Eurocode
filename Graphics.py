# Import Standard Libraries
import logging
# import scipy as np
import matplotlib.pyplot as plt

class TeeBeamPlot:
    bf = 1200
    hf = 150

    bw = 300
    h = 1200
    d_P = 957
    e_P = 43
    d_S = 1000

    Ap = 2400
    As = 0

    fck = 30
    fsk = 500
    fpk = 1860

    Ep = 195000
    Pm = 2772000
    units = "MPa"

    def __init__(self):
        self.sigma_p = self.Pm/self.Ap
        self.hw = self.h-self.hf
        Af = self.hf*self.bf
        Aw = self.bw*self.hw

        self.Ac = Af + Aw
        self.na = (Af*(self.hw+self.hf/2)+Aw*self.hw/2)/self.Ac
        self.I = self.bw*self.hw**3/12 + self.bf*self.hf**3/12 \
            + Af*(self.h-self.na-self.hf/2)**2 + Aw*(self.hw/2-self.na)
        logging.debug("    Ac = {:6.2f}, I = {:10.2f}".format(self.Ac, self.I))
        self.ta = self.h/100

    def plotTeeBeam(self):
        x0 = self.bf
        x1 = (self.bf + self.bw)/2
        x2 = (self.bf - self.bw)/2

        y0 = self.h
        y1 = self.h - self.hf
        X = [0, x0, x0, x1, x1, x2, x2, 0, 0]
        Y = [y0, y0, y1, y1, 0, 0, y1, y1, y0]
                
        self.fig = plt.figure(figsize=(9, 3))
        ax = self.fig.add_subplot(131)
        ax.plot(X, Y, linewidth=1, color='k')
        ax.fill_between(X, Y,alpha=0.5)

        # ax.plot([0, self.bf],[self.na, self.na], linestyle='--',color='k')

        ax.plot(self.bf/2, self.na, marker = '+', color='k')

        if self.Ap != 0:
            ax.plot(self.bf/2, self.h-self.d_P, marker = 'o', color='r')

        if self.As != 0:
            ax.plot(self.bf/2, self.h-self.d_S, marker = 'o', color='k')

        ax.plot(self.bf/2, self.na+self.e_P, marker='x', color='k')
        # plt.show()

    def plotStressDist(self,Xu,*args):

        ax = self.fig.add_subplot(132)
        ax.plot([0,0], [0,self.h], linewidth=2, color='k')

        # Moment Capacity
        if Xu < self.hf:
            Nc1 = args[0]
            print(Nc1)
        elif Xu > self.hf:


            [fcd, Nc1,y1,Nc2,y2,Ns,y3] = args

            h0 = self.h - self.hf
            self.__plotConcreteStress(ax,fcd,self.h,h0)

            h1 = self.h - Xu
            self.__plotConcreteStress(ax, fcd, self.h, h1)

            self.__plotSteelStress(ax, Ns, y3)
            print(Ns)
            # y1 = [self.h, self.h, h1, h1, self.h]


            h1 = h1 + Xu/2

            plt.xlim((-1.2*fcd,1.2*fcd))
            plt.show()

    # Private:

    def __plotConcreteStress(self,ax,fcd,h_top,h_bot):
        X0 = [-fcd, 0]
        Y0 = [h_top, h_top]
        Y1 = [h_bot, h_bot]
        ax.fill_between(X0, Y0, Y1, color='b', alpha=0.3)
        ax.arrow(fcd, (h_top+h_bot)/2, -9*fcd/10, 0,
        width=self.ta,  fc='b', ec='b', head_width=3*self.ta, head_length=self.ta/10)
    
    def __plotSteelStress(self,ax,Ns,h):
        if Ns > 0:
            ax.arrow(0, h, 9*Ns/10, 0,
            width=self.ta,  fc='r', ec='r', head_width=3*self.ta, head_length=self.ta/10)

if __name__ == "__main__":
    beam = TeeBeamPlot()
    beam.plotTeeBeam()
    Xu = 224
    fcd = 20
    Nc1 = 27
    y1 = 282
    Nc2 = 107
    y2 = 267
    Ns = 2
    y3 = 18

    beam.plotStressDist(Xu,fcd,Nc1,y1,Nc2,y2,Ns,y3)

