# Import Standard Liself.braries
import logging
import scipy as np
from utilities import *

# clself.Ass Concreteself.beam():

# self.Assume steel yields self.before concrete crushes


class RectangularBeam():
    h = 609.6
    d = 546.1
    b = 558.8
    As = 285.75
    fck = 27.6
    fy = 414
    rho = 25
    Es = 200000
    units = "MPa"
    # h = 24
    # d = 21.5
    # b = 22
    # As = 6
    # fck = 4000
    # fy = 60
    # rho = 145
    # Es = 29000000
    # units = "psi"

    def ACI_cracking_moment(self):
        """ Mcr = moment capacity prior to cracking """
        Ec = ACI_elastic_modulus(self.fck, self.rho, self.units)
        fr = ACI_tensile_strength(self.fck, self.units)
        n = self.Es/Ec
        # Area of uncracked concrete
        Ac = self.b*self.h
        # Equivalent concrete area of steel
        As = (n-1)*self.As
        # Distance from neutral axis to bottom of beam
        y_bott = (Ac*self.h/2 + As*(self.h-self.d))/(Ac+As)
        # Moment of inertia of uncracked concrete
        I_uncr = self.b*self.h**3/12\
            + Ac*(self.h/2-y_bott)**2\
            + As*(y_bott-(self.h-self.d))**2
        # Moment capacity prior to cracking
        Mcr = fr*I_uncr/y_bott
        return Mcr

    def ACI_design_moment(self):
        """ Mu = ultimate moment """
        beta1 = ACI_beta(self.fck, self.units)
        a = self.As*self.fy/(self.b*0.85*self.fck)
        Xu = a/beta1

        return self.As*self.fy*(self.d-a/2)

    def EC2_cracking_moment(self):
        Ec = EC2_elastic_modulus(self.fck, self.units)
        fr = EC2_flex_tensile_strength(self.fck,self.h,self.units)
        logging.debug("    fr = {} {}".format(fr,self.units))
        n = self.Es/Ec
        logging.debug("    modular ratio = {}".format(n))
        # Area of uncracked concrete
        Ac = self.b*self.h
        # Equivalent concrete area of steel
        As = (n-1)*self.As
        # Distance from neutral axis to bottom of beam
        y_bott = (Ac*self.h/2 + As*(self.h-self.d))/(Ac+As)
        logging.debug(y_bott)
        # Moment of inertia of uncracked concrete
        I_uncr = self.b*self.h**3/12\
            + Ac*(self.h/2-y_bott)**2\
            + As*(y_bott-(self.h-self.d))**2
        logging.debug(I_uncr)
        # Moment capacity prior to cracking
        Mcr = fr*I_uncr/y_bott
        return Mcr

    def EC2_design_moment(self):
        [alpha, self.beta] = EC2_alpha_beta(self.fck, self.units)
        fcd = self.fck/1.5
        ecu3 = EC2_max_strain(self.fck)
        Xu = self.As*self.fy/(alpha*self.b*fcd)
        return self.As*self.fy*(self.d-self.beta*Xu)

    def JSCE_moment_capacity(self):
        k1 = min(1 - 0.003*self.fck, 0.85)
        ecu = min((155-self.fck)/30000, 0.0035)
        self.beta = 0.52*80*ecu
        fcd = self.fck/1.5

        a = self.As*self.fy/(self.b*k1*fcd)
        Xu = a/self.beta

        return self.As*self.fy*(self.d-a/2)


if __name__ == "__main__":

    logging.basicConfig(level=logging.DEBUG, format='%(message)s')
    beam = RectangularBeam()

    MRd = beam.ACI_cracking_moment()
    logging.debug(MRd/12000)
    MRd = beam.EC2_cracking_moment()
    logging.debug(MRd/12000)
