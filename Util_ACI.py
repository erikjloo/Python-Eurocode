# Import Standard Libraries
import logging
import scipy as np

# Import Local Libraries
from Utilities import *

#===========================================================================
#   ACI Equations - Material properties
#===========================================================================


def elastic_modulus(fck, rho=145, units="psi"):
    """ Input:  fck = char. comp. strength of concrete
                rho = density of concrete (default = 145)
                units = "MPa" or "psi" (default = "psi") 
        Output: Ec = mean elastic modulus of concrete """
    fck = convert_2_psi(fck, units)
    rho = convert_2_pcf(rho, units)
    Ec = 33*rho**1.5*np.sqrt(fck)
    return Ec if units == "psi" else convert_2_MPa(Ec, "psi")


def tensile_strength(fck, units="psi"):
    """ Input:  fck = char. comp. strength of concrete 
                units = "MPa" or "psi" (default = "psi")
        Output: fctm = mean tensile strength of concrete """
    fck = convert_2_psi(fck, units)
    la = 1  # lambda factor necessary
    fctm = 7.5*la*np.sqrt(fck)
    return fctm if units == "psi" else convert_2_MPa(fctm, "psi")


def ultimate_strain(fck, units="psi"):
    """ Input:  fck = char. comp. strength of concrete
                units = "MPa" or "psi" (default = "psi") 
        Output: ecu = ultimate strain of concrete """
    fck = convert_2_psi(fck, units)
    return 3/1000

#===========================================================================
#   ACI Equations - Parameters
#===========================================================================


def beta(fck, units="psi"):
    """ Input:  fck = char. comp. strength of concrete
                units = "MPa" or "psi" (default = "psi") 
        Output: beta_1 = a/c in Whitney stress block """
    fck = convert_2_psi(fck, units)
    beta1 = 0.85-0.05*(fck-4000)/1000
    return min(max(beta1, 0.65), 0.85)


#===========================================================================
#   ACI Equations - Maximum reinforcement (Ductility)
#===========================================================================

def ductility_requirement(c, d, type="beam"):
    """ Input:  e_T = avg. tensile strain of steel
                type = "beam", "ties" or "spiral"
        Output: phi = strength reduction factor """
    e_T = (0.003/c)*(d-c)
    if type == "beam" and e_T >= 0.004:
        phi = 0.65 + 0.25/0.003*(e_T - 0.002)
        logging.info("    e_T = {:6.5f} > e_T,min = 0.004. OK".format(e_T))
        return min(phi, 0.9)
    elif type == "ties":
        phi = 0.65 + 0.25/0.003*(e_T - 0.002)
        return min(max(phi, 0.65), 0.9)
    elif type == "spiral":
        phi = 0.75 + 0.15/0.003*(e_T - 0.002)
        return min(max(phi, 0.75), 0.9)
    else:
        logging.info("    e_T = {:6.5f} < e_T,min = 0.004. Not OK".format(e_T))

#===========================================================================
#   ACI Equations - Minimum reinforcement (Md > Mcr)
#===========================================================================

def steel_ratio(As, fck, fyk, b, d, units="psi"):
    """ Input:  As = area of reinforcement steel
                fck = char. comp. strength of concrete
                fyk = char. yield stress of reinforcement
                b = width of beam portion in compression
                d = distance from comp. to reinforcement 
                units = "MPa" or "psi" (default = "psi")
        Output: A_min = minimum reinforcment area
                A_max = maximum reinforcement area """
    [fck, fyk] = convert_2_psi(np.array([fck, fyk]), units)
    [b, d] = convert_2_in(np.array([b, d]), units)
    As = convert_2_in2(As, units)

    A_min = max(3*np.sqrt(fck)/fyk, 200/fyk) * (b*d)
    A_max = 0.365*(fck/fyk)*beta(fck) * (b*d)

    compare_steel_area(As, A_min, A_max)

    return [A_min, A_max] if units == "psi" else [A_min*25.4**2, A_max*25.4**2]

#===========================================================================
#   ACI Equations - Design Aids
#===========================================================================


def barsize(Nbar, units="psi"):
    """ bar_area = barsize(Nbar or bar_area) """
    bar_area = {3: 0.11, 4: 0.20, 5: 0.31, 6: 0.44, 7: 0.60,
                8: 0.79, 9: 1, 10: 1.27, 11: 1.56, 14: 2.25, 18: 4.00}
    return bar_area.get(Nbar, Nbar)


def min_beam_width(numbar, rebar, units="psi"):
    """ min_width = min_beam_width(numbar,rebar) """
    if numbar > 1 & numbar < 7 & units == "psi":
        BW = [[7.0, 8.5, 10.0, 11.5, 13.0],
              [7.0, 8.5, 10.5, 12.0, 13.5],
              [7.0, 9.0, 11.0, 12.5, 14.0],
              [7.5, 9.0, 11.0, 13.0, 15.0],
              [7.5, 9.5, 11.5, 13.5, 15.5],
              [8.0, 10.0, 12.5, 14.5, 17.0],
              [8.0, 10.5, 13.0, 15.5, 18.0],
              [8.5, 11.0, 14.0, 17.0, 19.5]]
        return BW[rebar-4, numbar-2]
    else:
        raise ValueError()


