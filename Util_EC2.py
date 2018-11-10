# Import Standard Libraries
import logging
import scipy as np

# Import Local Libraries
from Utilities import *

#===========================================================================
#   EC2 Equations - Material properties
#===========================================================================


def elastic_modulus(fck, units="MPa"):
    """ Input:  fck = char. comp. strength of concrete
                units = "MPa" or "psi" (default = "MPa")
        Output: Ec = mean elastic modulus of concrete """
    fck = convert_2_MPa(fck, units)
    fcm = fck+8
    Ec = 22000*(fcm/10)**0.3
    return Ec if units == "MPa" else convert_2_psi(Ec, "MPa")


def tensile_strength(fck, units="MPa"):
    """ Input:  fck = char. comp. strength of concrete
                units = "MPa" or "psi" (default = "MPa")
        Output: fctm = mean tensile strength of concrete """
    fck = convert_2_MPa(fck, units)
    fcm = fck+8
    fctm = 0.3*fck**(2/3) if fck <= 50 else 2.12*np.log(1+fcm/10)
    return fctm if units == "MPa" else convert_2_psi(fctm, "MPa")


def flex_tensile_strength(fck, h, units="MPa"):
    """ Input:  fck = char. comp. strength of concrete
                h = height of reinforced concrete beam
                units = "MPa" or "psi" (default = "MPa") 
        Output: fctm,fl = mean tensile strength for flexure """
    fck = convert_2_MPa(fck, units)
    fctm = tensile_strength(fck)
    h = convert_2_mm(h, units)
    fctm = min((1.6-h/1000)*fctm, fctm)
    return fctm if units == "MPa" else convert_2_psi(fctm, "MPa")


def ultimate_strain(fck, units="MPa"):
    """ Input:  fck = char. comp. strength of concrete
                units = "MPa" or "psi" (default = "MPa") 
        Output: ecu3 = ultimate tensile strain """
    fck = convert_2_MPa(fck, units)
    ecu3 = 2.6+35*((90-fck)/100)**4
    return min(ecu3, 3.5)/1000

#===========================================================================
#   EC2 Equations - Parameters
#===========================================================================


def alpha_beta(fck, units="MPa"):
    """ Input:  fck = char. comp. strength of concrete
                units = "MPa" or "psi" (default = "MPa") 
        Output: alpha = factor for bilinear stress block
                beta = (dist. from comp. to Nc)/Xu """
    fck = convert_2_MPa(fck, units)
    alpha = np.ceil((9E-05*fck**2 - 0.0177*fck + 1.4032)*100)/100
    beta = np.ceil((4E-05*fck**2 - 0.0071*fck + 0.634)*100)/100
    return [min(alpha, 0.75), min(beta, 0.39)]


def lambda_eta(fck, units="MPa"):
    """ Input:  fck = char. comp. strength of concrete
                units = "MPa" or "psi" (default = "MPa") 
        Output: la = (height of compressive zone)/Xu
                eta = factor for "Whitney" stress block """
    fck = convert_2_MPa(fck, units)
    la = min(0.8-(fck-50)/400, 0.8)
    eta = min(1-(fck-50)/200, 1.0)
    return [la, eta]

#===========================================================================
#   EC2 Equations - Maximum reinforcement (Ductility)
#===========================================================================


def ductility_requirement(Xu, d, fck, fyd, units="MPa"):
    """ Input:  Xu = dist. from comp. to neutral axis 
                d = dist. from comp. to reinforcement
                fck = char. comp. strength of concrete
                fyd = design steel yield stress
                units = "MPa" or "psi" (default = "MPa")
        Output: Xu_max = Max. dist. to neutral axis """
    [fck, fyd] = convert_2_MPa(np.array([fck, fyd]), units)
    
    ecu = ultimate_strain(fck)  # units="MPa"
    Xu_max = min(ecu*10**6/(ecu*10**6+7*fyd), 0.535)*d
    if Xu < Xu_max:
        logging.info(
            "    Xu = {:6.2f} < Xu_max = {:6.2f}. OK".format(Xu, Xu_max))
    else:
        logging.info(
            "    Xu = {:6.2f} > Xu_max = {:6.2f}. Not OK".format(Xu, Xu_max))
    return Xu_max

#===========================================================================
#   EC2 Equations - Minimum reinforcement (Md > Mcr)
#===========================================================================


def steel_ratio(As, fck, fyk, b, d, h, Xu, units="MPa"):
    """ Input:  As = area of reinforcement steel
                fck = char. comp. strength of concrete
                fyk = char. yield stress of reinforcement
                b = width of beam portion in compression
                d = dist. from comp. to reinforcement
                h = height of reinforced concrete beam
                Xu = maximum dist. to neutral axis
                units = "MPa" or "psi" (default = "MPa")
        Output: A_min = minimum reinforcement area
                A_max = maximum reinforcement area """
    [fck, fyk] = convert_2_MPa(np.array([fck, fyk]), units)
    [b, d, h, Xu] = convert_2_mm(np.array([b, d, h, Xu]), units)
    As = convert_2_mm2(As, units)

    fctm = flex_tensile_strength(fck, h)  # units="MPa"
    A_min = max((0.26*fctm/fyk),0.0013)* (b*d)
    

    fcd = fck/1.5
    fyd = fyk/1.15

    alpha = alpha_beta(fck)[0]  # units="MPa"
    A_max = min(alpha*(fcd/fyd)*b*Xu, 0.4*b*d)

    compare_steel_area(As, A_min, A_max)

    return [A_min, A_max] if units == "MPa" else [A_min/(25.4**2), A_max/(25.4**2)]
