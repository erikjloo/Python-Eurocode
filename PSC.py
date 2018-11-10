# Import Standard Libraries
import logging
import scipy as np
import matplotlib.pyplot as plt

# Import Local Libraries
import Util_ACI as ACI
import Util_EC2 as EC2
from RC import ReinforcedBeam


#===========================================================================
#   Prestressed Beam
#===========================================================================

class PrestressedBeam(ReinforcedBeam):
    """ Abstract Reinforced Beam Class
    
    Static Members:
        gamma_c = 1.5
        gamma_s = 1.15
        gamma_P = 1.1
        niter = 20

    Instance Members:
    
    Public Methods:
    
    Private Methods:
        __plotArrow(self, ax, x, y, dx, dy, color)
        __plotConcreteStress(self, ax, fcd, h_top, h_bot) 
    """
    gamma_P = 1.1
    
    def __init__(self):
        ReinforcedBeam.__init__(self)

    def rotation_capacity(self):
        pass

#===========================================================================
#   Rectangular Beam
#===========================================================================


class RectangularBeam(PrestressedBeam):
    h = 1000
    d_P = 800
    e_P = 300
    b = 800
    Ap = 4300
    fck = 45
    fpk = 1860

    Ep = 195000
    Pm = 4800000
    units = "MPa"

    def __init__(self):

        self.sigma_p = self.Pm/self.Ap
        self.Ac = self.b*self.h
        self.I = self.b*self.h**3/12
        self.W = self.b*self.h**2/6

    def EC2_design_moment(self):
        [alpha, beta] = EC2.alpha_beta(self.fck, self.units)
        logging.debug("    a = {:3.2f}, b = {:3.2f}".format(alpha, beta))
        fcd = self.fck/self.gamma_c
        fp01d = 0.9*self.fpk/self.gamma_P
        fpd = self.fpk/self.gamma_P
        epy = fp01d/self.Ep
        logging.debug("    fcd = {:6.2f}, fpd = {:6.2f}".format(fcd, fpd))
        logging.debug("    fp01d = {:6.2f}, epy = {:6.5f}".format(fp01d, epy))

        ecu3 = EC2.ultimate_strain(self.fck, self.units)
        epu = 0.035

        # Assume fp = 0.95*fpd
        fp = 0.95*fpd
        logging.debug("    Assume fp = 0.95fpd = {:6.2f}".format(fp))
        for iter in range(self.niter):

            # Find Xu using horizontal equilibrium
            Xu = (self.Ap*(fp-self.sigma_p)+self.Pm)/(alpha*self.b*fcd)
            logging.debug(
                "Iteration {:1.0f}: \n    Xu = {:6.2f}".format(iter, Xu))

            # Is the assumption correct?
            Delta_ep = (self.d_P-Xu)*ecu3/Xu
            ep = self.sigma_p/self.Ep + Delta_ep
            fp_new = fp01d + (ep-epy)*(fpd-fp01d)/(epu-epy)
            logging.debug("    ep = {:6.5f}, fp = {:6.2f}".format(ep, fp_new))

            # Approximation Error
            if abs((fp_new-fp)/fp_new) < 0.01:
                break
            # New assumption fp = fp_new
            fp = fp_new

        # Moment Capacity
        Nc = alpha*self.b*Xu*fcd
        y1 = self.d_P - self.e_P - beta*Xu
        Delta_P = self.Ap*fp_new - self.Pm
        MRd = Nc*y1 + Delta_P*self.e_P
        logging.info("    MRd = {:5.2f}".format(MRd))
        return MRd

#===========================================================================
#   T Beam
#===========================================================================


class TeeBeam(PrestressedBeam):
    bf = 1200
    hf = 150

    bw = 300
    h = 1200

    d_S = 1000
    d_P = 1000
    d_P0 = 500

    As = 0
    Ap = 2400

    fck = 30
    fyk = 500
    fpk = 1860

    Es = 205000
    k = 0.01
    Ep = 195000
    Pm = 2772000
    units = "MPa"

    def __init__(self):

        PrestressedBeam.__init__(self)

        self.sigma_p = self.Pm/self.Ap
        self.f = self.d_P - self.d_P0
        hw = self.h-self.hf
        Af = self.hf*self.bf
        Aw = self.bw*hw

        self.Ac = Af + Aw
        self.yb = (Af*(hw+self.hf/2)+Aw*hw/2)/self.Ac
        self.yt = self.h - self.yb
        logging.debug("    yt = {:6.2f}, yb = {:}".format(self.yb,self.yt))

        self.e_P = self.d_P - self.yt
        print(self.e_P)
        
        self.I = self.bw*hw**3/12 + self.bf*self.hf**3/12 \
            + Af*(self.h-self.yb-self.hf/2)**2 + Aw*(hw/2-self.yb)
        logging.debug("    Ac = {:6.2f}, I = {:10.2f}".format(self.Ac, self.I))

        self.Wt = self.I/self.yt
        self.Wb = self.I/self.yb

        logging.debug("    Wt = {:6.2f}, Wb = {:6.2f}".format(self.Wt, self.Wb))
        

    def EC2_design_moment(self):
        [alpha, beta] = EC2.alpha_beta(self.fck, self.units)
        [la, eta] = EC2.lambda_eta(self.fck, self.units)
        logging.debug("    a = {:3.2f}, b = {:3.2f}".format(alpha, beta))
        logging.debug("    lambda = {:3.2f}, eta = {:3.2f}".format(la, eta))
        fcd = self.fck/self.gamma_c
        fyd = self.fyk/self.gamma_s
        fp01d = 0.9*self.fpk/self.gamma_P
        fpd = self.fpk/self.gamma_P

        epy = fp01d/self.Ep
        ey = fyd/self.Es

        logging.debug("    fcd = {:6.2f}, fpd = {:6.2f}".format(fcd, fpd))
        logging.debug("    fp01d = {:6.2f}, epy = {:6.5f}".format(fp01d, epy))

        ecu3 = EC2.ultimate_strain(self.fck, self.units)
        epu = 0.035

        # Assume fp = 0.95*fpd and fs = fyd
        fp = 0.95*fpd
        fs = fyd
        logging.debug("    Assume fp = 0.95fpd = {:6.2f}".format(fp))

        Xu = 300

        for iter in range(self.niter):

            if Xu < self.hf:
                Xu = (self.Ap*(fp-self.sigma_p)+self.Pm)/(alpha*self.bf*fcd)

            elif Xu > self.hf:
                X = (self.Ap*fp-eta*fcd*self.hf *
                     (self.bf-self.bw))/(eta*fcd*self.bw)
                Xu = X/la
            
            logging.debug("Iteration {:1.0f}: \n    Xu = {:6.2f}".format(iter,Xu))

            # Is the assumption correct?
            Delta_ep = (self.d_P-Xu)*ecu3/Xu
            ep = self.sigma_p/self.Ep + Delta_ep
            fp_new = fp01d + (ep-epy)*(fpd-fp01d)/(epu-epy)
            logging.debug("    ep = {:6.5f}, fp = {:6.2f}".format(ep, fp_new))

            es = (self.d_S-Xu)*ecu3/Xu
            fs_new = fyd + (es-ey)*self.k

            # Approximation Error
            
            if abs((fp_new-fp)/fp_new) < 0.001 and abs((fs_new-fs)/fs_new) < 0.001:
                break

            # New assumption f = f_new
            fp = fp_new
            fs = fs_new

        # Moment Capacity
        if Xu < self.hf:
            # _________
            #   |  |     ||
            #   |  |     ||   b*Xu
            #   |  |_    ||<------ Nc
            #   |    \_  ||    |
            #   |      \_||    y1
            #   |        ||    |
            #   dp ----- ||<------ Pm　＠　neutral axis
            #   |    |   ||    y2
            #   |    e   ||<------ Ns
            #   |    |   ||
            #  --------- ||<------Delta Pm
            #            ||

            Nc = alpha*self.bf*Xu*fcd
            y1 = self.d_P - self.e_P - beta*Xu
            logging.info("    Nc = {:6.2f}, y1 = {:6.2f}".format(Nc,y1))

            Ns = self.As*fs
            y2 = self.d_S - y1 - beta*Xu
            logging.info("    Ns = {:6.2f}, y2 = {:6.2f}".format(Ns,y2))
            
            Delta_P = self.Ap*fp - self.Pm
            logging.info(
                "    Delta_P = {:6.2f}, e_P = {:6.2f}".format(Delta_P, self.e_P))

            MRd = Nc*y1 + Ns*y2 + Delta_P*self.e_P
        
        elif Xu > self.hf:
            # _________
            #   |  |     ||
            #   |  |     ||<----------- Nc1
            #   |  |_____||<---- Nc2 |
            #   |  |     ||   |      |
            #   |  |_____||   y2     y1
            #   |        ||   |      |
            #   dp ----  ||<------------ Pm　＠　neutral axis
            #   |    |   ||    y3
            #   |    e   ||<------ Ns
            #   |    |   ||
            #  --------- ||<------Delta Pm
            #            ||

            Nc1 = (self.bf-self.bw)*self.hf*fcd
            y1 = self.d_P - self.e_P - self.hf/2
            logging.info("    Nc1 = {:6.2f}, y1 = {:6.2f}".format(Nc1, y1))

            Nc2 = self.bw*X*fcd
            y2 = self.d_P - self.e_P - X/2
            logging.info("    Nc2 = {:6.2f}, y2 = {:6.2f}".format(Nc2, y2))

            Ns = self.As*fs
            y3 = self.d_P - y1 - self.hf/2
            logging.info("    Ns = {:6.2f}, y3 = {:6.2f}".format(Ns, y3))

            Delta_P = self.Ap*fp - self.Pm
            logging.info("    Delta_P = {:6.2f}, e_P = {:6.2f}".format(Delta_P, self.e_P))

            MRd = Nc1*y1 + Nc2*y2 + Ns*y3 + Delta_P*self.e_P
        logging.info("    MRd = {:5.2f}".format(MRd))
        return MRd

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
        ax.fill_between(X, Y, alpha=0.5)

        ax.plot([0, self.bf],[self.yb, self.yb], linestyle='--',color='k')

        ax.plot(self.bf/2, self.yb, marker='+', color='k')

        if self.Ap != 0:
            ax.plot(self.bf/2, self.h-self.d_P, marker='o', color='r')

        if self.As != 0:
            ax.plot(self.bf/2, self.h-self.d_S, marker='o', color='k')

        ax.plot(self.bf/2, self.h-self.d_P0, marker='o', color='r')
        plt.show()

if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG, format='%(message)s')
    beam = TeeBeam()
    beam.EC2_design_moment()
    beam.plotTeeBeam()
