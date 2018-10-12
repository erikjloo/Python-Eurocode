# Import Standard Liself.braries
import logging
import scipy as np
from Util_ACI import *
from Util_EC2 import *

# clself.Ass Concreteself.beam():

# self.Assume steel yields self.before concrete crushes


# class PrestressedBeam():

#===========================================================================
#   Rectangular Beam
#===========================================================================

class RectangularBeam():
    # h = 609.6
    # d = 546.1
    # b = 558.8
    # As = 285.75
    # fck = 27.6
    # fyk = 414
    # rho = 25
    # Es = 200000
    # units = "MPa"
    h = 24
    d = 21.5
    b = 22
    As = 6
    fck = 4000
    fyk = 60000
    rho = 145
    Es = 29000000
    units = "psi"

    #---------------------------------------------------------------------------
    #   ACI Equations
    #---------------------------------------------------------------------------

    def ACI_cracking_moment(self):
        """ Mcr = moment capacity prior to cracking """
        logging.debug("Uncracked moment capacity per ACI")
        Ec = ACI_elastic_modulus(self.fck, self.rho, self.units)
        fr = ACI_tensile_strength(self.fck, self.units)
        n = self.Es/Ec
        logging.debug("    fr = {:6.2f}, Es/Ec = {:3.2f}".format(fr, n))
        Ac = self.b*self.h
        logging.debug("    Area of uncracked concrete = {:6.2f}".format(Ac))
        As = (n-1)*self.As
        logging.debug("    Equivalent area of steel = {:6.2f}".format(As))
        y_bott = (Ac*self.h/2 + As*(self.h-self.d))/(Ac+As)
        logging.debug(
            "    Dist. NA to bottom: y_bott = {:6.2f}".format(y_bott))
        I_uncr = self.b*self.h**3/12\
            + Ac*(self.h/2-y_bott)**2\
            + As*(y_bott-(self.h-self.d))**2
        logging.debug("    Moment of inertia: = {:10.2f} ".format(I_uncr))
        return fr*I_uncr/y_bott

    def ACI_elastic_moment(self):
        pass

    def ACI_design_moment(self):
        """ Mu = ultimate moment """
        beta1 = ACI_beta(self.fck, self.units)
        logging.debug("    beta1 = {:3.2f}".format(beta1))
        c = (self.As*self.fyk)/(0.85*self.fck*self.b)/beta1
        logging.debug("    c = {:6.2f}".format(c))
        phi = ACI_ductility_requirement(c, self.d, type="beam")
        logging.debug("    phi = {:3.2f}".format(phi))
        ACI_steel_ratio(self.As, self.fck, self.fyk,
                        self.b, self.d, self.units)
        return phi*self.As*self.fyk*(self.d-beta1*c/2)

    #---------------------------------------------------------------------------
    #   EC2 Equations
    #---------------------------------------------------------------------------

    def EC2_cracking_moment(self):
        logging.debug("Uncracked moment capacity per EC2")
        Ec = EC2_elastic_modulus(self.fck, self.units)
        fr = EC2_flex_tensile_strength(self.fck, self.h, self.units)
        n = self.Es/Ec
        logging.debug("    fr = {:6.2f}, Es/Ec = {:3.2f}".format(fr, n))
        Ac = self.b*self.h
        logging.debug("    Area of uncracked concrete = {:6.2f}".format(Ac))
        As = (n-1)*self.As
        logging.debug("    Equivalent area of steel = {:6.2f}".format(As))
        y_bott = (Ac*self.h/2 + As*(self.h-self.d))/(Ac+As)
        logging.debug(
            "    Dist. NA to bottom: y_bott = {:6.2f}".format(y_bott))
        I_uncr = self.b*self.h**3/12\
            + Ac*(self.h/2-y_bott)**2\
            + As*(y_bott-(self.h-self.d))**2
        logging.debug("    Moment of inertia: = {:10.2f} ".format(I_uncr))
        return fr*I_uncr/y_bott

    def EC2_design_moment(self):
        [alpha, beta] = EC2_alpha_beta(self.fck, self.units)
        logging.debug(
            "    alpha = {:3.2f}, beta = {:3.2f}".format(alpha, beta))
        fyd = self.fyk/1.15
        fcd = self.fck/1.5
        logging.debug("    fyd = {:6.2f}, fcd = {:6.2f}".format(fyd, fcd))
        Xu = self.As*fyd/(alpha*self.b*fcd)
        Xu_max = EC2_ductility_requirement(
            Xu, self.d, self.fck, fyd, self.units)

        EC2_steel_ratio(self.As, self.fck, self.fyk, self.b,
                        self.d, self.h, Xu_max, self.units)

        return self.As*fyd*(self.d-beta*Xu)

    def JSCE_moment_capacity(self):
        k1 = min(1 - 0.003*self.fck, 0.85)
        ecu = min((155-self.fck)/30000, 0.0035)
        beta = 0.52*80*ecu
        fcd = self.fck/1.5

        a = self.As*self.fyk/(self.b*k1*fcd)
        Xu = a/beta

        return self.As*self.fyk*(self.d-a/2)


class DoublyReinforcedBeam():
    a = 2


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG, format='%(message)s')

    print("\n Ec for C30/37: \n")
    print(ACI_elastic_modulus(30, rho=25, units="MPa"))
    print(EC2_elastic_modulus(30))
    print("\n Tensile strength for C30/37: \n")
    print(ACI_tensile_strength(30, units="MPa"))
    print(EC2_tensile_strength(30))
    print("\n Tensile strength for C55/67: \n")
    print(ACI_tensile_strength(55, units="MPa"))
    print(EC2_tensile_strength(55))
    print("\n Ultimate strain for C30/37 \n")
    print(ACI_ultimate_strain(30, units="MPa"))
    print(EC2_ultimate_strain(30))
    print("\n Ultimate strain for C55/67 \n")
    print(ACI_ultimate_strain(55, units="MPa"))
    print(EC2_ultimate_strain(55))
    print("\n alpha & Beta factors: \n")
    print(EC2_alpha_beta(10))
    print(EC2_alpha_beta(60))
    print(EC2_alpha_beta(80))
    print("\n Lambda & Eta factors: \n")
    print(EC2_lambda_eta(10))
    print(EC2_lambda_eta(60))
    print(EC2_lambda_eta(80))

    beam = RectangularBeam()
    print('\n Cracking Moment: \n')
    MRd = beam.ACI_cracking_moment()
    print(MRd/12000)
    MRd = beam.EC2_cracking_moment()
    print(MRd/12000)

    print('\n Design Moment: \n')
    MRd = beam.ACI_design_moment()
    print(MRd/12000)
    MRd = beam.EC2_design_moment()
    print(MRd/12000)
