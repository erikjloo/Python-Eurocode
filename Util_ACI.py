# Import Standard Libraries
import logging
import scipy as np

# Import Local Libraries
from Utilities import *

#===========================================================================
#   ACI Equations - Material properties
#===========================================================================

def ACI_elastic_modulus(fck, rho=145, units="psi"):
    """ Input:  fck = char. comp. strength of concrete
                rho = density of concrete (default = 145)
                units = "MPa" or "psi" (default = "psi") 
        Output: Ec = mean elastic modulus of concrete """
    fck = convert_2_psi(fck, units)
    rho = convert_2_pcf(rho, units)
    Ec = 33*rho**1.5*np.sqrt(fck)
    return Ec if units == "psi" else convert_2_MPa(Ec, "psi")


def ACI_tensile_strength(fck, units="psi"):
    """ Input:  fck = char. comp. strength of concrete 
                units = "MPa" or "psi" (default = "psi")
        Output: fctm = mean tensile strength of concrete """
    fck = convert_2_psi(fck, units)
    la = 1  # lambda factor necessary
    fctm = 7.5*la*np.sqrt(fck)
    return fctm if units == "psi" else convert_2_MPa(fctm, "psi")


def ACI_ultimate_strain(fck, units="psi"):
    """ Input:  fck = char. comp. strength of concrete
                units = "MPa" or "psi" (default = "psi") 
        Output: ecu = ultimate strain of concrete """
    fck = convert_2_psi(fck, units)
    return 3/1000

#===========================================================================
#   ACI Equations - Parameters
#===========================================================================


def ACI_beta(fck, units="psi"):
    """ Input:  fck = char. comp. strength of concrete
                units = "MPa" or "psi" (default = "psi") 
        Output: beta_1 = a/c in Whitney stress block """
    fck = convert_2_psi(fck, units)
    beta1 = 0.85-0.05*(fck-4000)/1000
    return min(max(beta1, 0.65), 0.85)


#===========================================================================
#   ACI Equations - Minimum and maximum reinforcement
#===========================================================================

def ACI_ductility_requirement(c, d, type="beam"):
    """ Input:  e_T = avg. tensile strain of steel
                type = "beam", "ties" or "spiral"
        Output: phi = strength reduction factor """
    e_T = (0.003/c)*(d-c)
    if type == "beam" and e_T >= 0.004:
        phi = 0.65 + 0.25/0.003*(e_T - 0.002)
        logging.info("    e_T = {:8.6f} > 0.004. OK".format(e_T))
        return min(phi, 0.9)
    elif type == "ties":
        phi = 0.65 + 0.25/0.003*(e_T - 0.002)
        return min(max(phi, 0.65), 0.9)
    elif type == "spiral":
        phi = 0.75 + 0.15/0.003*(e_T - 0.002)
        return min(max(phi, 0.75), 0.9)
    else:
        logging.error("    e_T = {:8.6f} < 0.004. Not OK".format(e_T))


def ACI_steel_ratio(As, fck, fyk, b, d, units="psi"):
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
    A_max = 0.365*(fck/fyk)*ACI_beta(fck) * (b*d)

    if As > A_min:
        logging.info(
            "    As = {:8.6f} > A_min = {:8.6f}. Md > Mcr".format(As, A_min))
    else:
        logging.error(
            "    As = {:8.6f} < A_min = {:8.6f}. Md < Mcr".format(As, A_min))

    return [A_min, A_max] if units == "psi" else [A_min, A_max]*(25.4**2)


def ACI_barsize(Nbar, units="psi"):
    """ bar_area = ACI_barsize(Nbar or bar_area) """
    bar_area = {3: 0.11, 4: 0.20, 5: 0.31, 6: 0.44, 7: 0.60,
                8: 0.79, 9: 1, 10: 1.27, 11: 1.56, 14: 2.25, 18: 4.00}
    return bar_area.get(Nbar, Nbar)


def ACI_min_beam_width(numbar, rebar, units="psi"):
    """ min_width = ACI_min_beam_width(numbar,rebar) """
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

#===========================================================================
#   EC2 Equations - Material properties
#===========================================================================


def EC2_elastic_modulus(fck, units="MPa"):
    """ Input:  fck = char. comp. strength of concrete
                units = "MPa" or "psi" (default = "MPa")
        Output: Ec = elastic modulus of concrete """
    fcm = convert_2_MPa(fck, units)
    fcm = fck+8
    Ec = 22000*(fcm/10)**0.3
    return Ec if units == "MPa" else convert_2_psi(Ec, "MPa")


def EC2_tensile_strength(fck, units="MPa"):
    """ Input:  fck = char. comp. strength of concrete
                units = "MPa" or "psi" (default = "MPa")
        Output: fctm = mean tensile strength of concrete """
    fck = convert_2_MPa(fck, units)
    fcm = fck+8
    fctm = 0.3*fck**(2/3) if fck <= 50 else 2.12*np.log(1+fcm/10)
    return fctm if units == "MPa" else convert_2_psi(fctm, "MPa")


def EC2_flex_tensile_strength(fck, h, units="MPa"):
    """ Input:  fck = char. comp. strength of concrete
                h = height of reinforced concrete beam
                units = "MPa" or "psi" (default = "MPa") 
        Output: fctm,fl = mean tensile strength for flexure """
    fctm = EC2_tensile_strength(fck, units)
    h = convert_2_mm(h, units)
    return min((1.6-h/1000)*fctm, fctm)


def EC2_ultimate_strain(fck, units="MPa"):
    """ Input:  fck = char. comp. strength of concrete
                units = "MPa" or "psi" (default = "MPa") 
        Output: ecu3 = ultimate tensile strain """
    fck = convert_2_MPa(fck, units)
    ecu3 = 2.6+35*((90-fck)/100)**4
    return min(ecu3, 3.5)/1000

#===========================================================================
#   EC2 Equations - Parameters
#===========================================================================


def EC2_alpha_beta(fck, units="MPa"):
    """ [alpha,beta] = EC_alpha_beta(fc) """
    fck = convert_2_MPa(fck, units)
    alpha = np.ceil((9E-05*fck**2 - 0.0177*fck + 1.4032)*100)/100
    beta = np.ceil((4E-05*fck**2 - 0.0071*fck + 0.634)*100)/100
    return [min(alpha, 0.75), min(beta, 0.38)]


def EC2_lambda_eta(fck, units="MPa"):
    """ [lambda, eta] = EC2_lambda_eta(fck) """
    fck = convert_2_MPa(fck, units)
    if fck > 50:
        la = 0.8-(fck-50)/400
        eta = 1-(fck-50)/200
    else:
        la = 0.8
        eta = 1.0
    # la = 0.8-(fck-50)/400
    # eta = 1-(fck-50)/200
    # la = min(la,0.8), min(eta, )
    return [la, eta]

#===========================================================================
#   EC2 Equations - Minimum and maximum reinforcement
#===========================================================================


def EC2_ductility_requirement(Xu, d, fck, fyd, units="MPa"):
    """ Input:  Xu = Distance from comp. to neutral axis 
                d = Distance from comp. to reinforcement
                fck = char. comp. strength of concrete
                fyd = design steel yield stress
                units = "MPa" or "psi" (default = "MPa")
        Output: Xu_max = Max. dist. to neutral axis """
    [fck, fyd] = convert_2_MPa(np.array([fck,fyd]),units)
    ecu = EC2_ultimate_strain(fck, units)
    Xu_max = min(ecu*10**6/(ecu*10**6+7*fyd), 0.535)*d
    if Xu < Xu_max:
        logging.info(
            "    Xu = {:8.6f} < Xu_max = {:8.6f}. OK".format(Xu, Xu_max))
    else:
        logging.error(
            "    Xu = {:8.6f} > Xu_max = {:8.6f}. Not OK".format(Xu, Xu_max))
    return Xu_max


def EC2_steel_ratio(As, fck, fyk, b, d, h, Xu_max, units="MPa"):
    """ Input:  As = area of reinforcement steel
                fck = char. comp. strength of concrete
                fyk = char. yield stress of reinforcement
                b = width of beam portion in compression
                d = distance from comp. to reinforcement
                h = height of reinforced concrete beam
                units = "MPa" or "psi" (default = "MPa")
        Output: A_min = minimum reinforcment area
                A_max = maximum reinforcement area """
    [fck, fyk] = convert_2_MPa(np.array([fck, fyk]), units)
    [b, d, h] = convert_2_mm(np.array([b, d, h]), units)
    As = convert_2_mm2(As, units)

    fctm = EC2_flex_tensile_strength(fck, h, units)
    rho_min = 0.26*fctm/fyk
    A_min = rho_min*b*d

    fcd = fck/1.5
    fyd = fyk/1.15

    alpha = EC2_alpha_beta(fck, "MPa")[0]
    print(alpha)
    A_max = alpha*(fcd/fyd)*b*Xu_max
    A_max = min(A_max, 0.4*b*d)

    if As > A_min:
        logging.info(
            "    As = {:8.6f} > A_min = {:8.6f}. Md > Mcr".format(As, A_min))
    else:
        logging.error(
            "    As = {:8.6f} < A_min = {:8.6f}. Md < Mcr".format(As, A_min))

    if As < A_max:
        logging.info(
            "    As = {:8.6f} > A_min = {:8.6f}. Md > Mcr".format(As, A_min))
    else:
        logging.error(
            "    As = {:8.6f} < A_min = {:8.6f}. Md < Mcr".format(As, A_min))

    return [A_min, A_max] if units == "MPa" else [A_min/(25.4**2), A_max/(25.4**2)]




if __name__ == "__main__":
    print("\n Ec for C30/37: \n")
    print(ACI_elastic_modulus(30, units="MPa"))
    print(EC2_elastic_modulus(30))
    print("\n Tensile strength for C30/37: \n")
    print(ACI_tensile_strength(30, units="MPa"))
    print(EC2_tensile_strength(30))
    print("\n Tensile strength for C55/67: \n")
    print(ACI_tensile_strength(55, units="MPa"))
    print(EC2_tensile_strength(55))
    print("\n Ultimate strain for C30/37 \n")
    print(EC2_ultimate_strain(30))
    print("\n Ultimate strain for C55/67 \n")
    print(EC2_ultimate_strain(55))
    print("\n Alpha & Beta factors: \n")
    print(EC2_alpha_beta(10))
    print(EC2_alpha_beta(60))
    print(EC2_alpha_beta(80))
    print("\n Lambda & Eta factors: \n")
    print(EC2_lambda_eta(10))
    print(EC2_lambda_eta(60))
    print(EC2_lambda_eta(80))
