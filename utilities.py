# Import Standard Libraries
import logging
import scipy as np


#===========================================================================
#   ACI Equations - Material Properties
#===========================================================================

def ACI_elastic_modulus(fck, rho=145, units="psi"):
    """ Input:  fck = char. comp. strength of concrete
                rho = density of concrete 
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
    fctm = 7.5*np.sqrt(fck)
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
    if fck < 2500:
        raise ValueError(
            "fck = {} {} is too small".format(fck, units))
    beta1 = 0.85-0.05*(fck-4000)/1000
    return min(max(beta1, 0.65), 0.85)


def ACI_steel_ratio(fck, fy, units="psi"):
    """ [ratio_min, ratio_max] = ACI_steel_ratio(fc, fy, units="psi")"""
    fck = convert_2_psi(fck, units)
    fy = convert_2_psi(fy, units)

    a1 = 3*np.sqrt(fck)/fy
    a2 = 200/fy
    ratio_min = max(a1, a2)
    beta1 = ACI_beta(fck)
    ratio_max = 0.365*fck*beta1/fy
    return [ratio_min, ratio_max]


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
#   EC2 Equations - Material Properties
#===========================================================================


def EC2_elastic_modulus(fck, units="MPa"):
    """ Input:  fck = char. comp. strength of concrete
                units = "MPa" or "psi" (default = "MPa")
        Output: Ec = elastic modulus of concrete """
    fck = convert_2_MPa(fck, units)
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


def EC2_flex_tensile_strength(fck,h,units="MPa"):
    """ Input:  fck = char. comp. strength of concrete
                units = "MPa" or "psi" (default = "psi") 
        Output: fctm,fl = tensile strength for flexure """
    fctm = EC2_tensile_strength(fck, units)
    h = convert_2_mm(h, units)
    return min((1.6-h/1000)*fctm,fctm)

def EC2_ultimate_strain(fck, units="MPa"):
    """ Input:  fck = char. comp. strength of concrete
                units = "MPa" or "psi" (default = "psi") 
        Output: ecu3 = ultimate tensile strain """
    fck = convert_2_MPa(fck, units)
    ecu3 = 2.6+35*((90-fck)/100)**4
    return min(ecu3,3.5)/1000

#===========================================================================
#   EC2 Equations - Parameters
#===========================================================================

def EC2_alpha_beta(fck, units="MPa"):
    """ [alpha,beta] = EC_alpha_beta(fc) """
    fck = convert_2_MPa(fck, units)
    if fck > 90:
        raise ValueError('fck = {} {} > 90 Mpa'.format(fck, units))
    alpha = np.ceil((9E-05*fck**2 - 0.0177*fck + 1.4032)*100)/100
    beta = np.ceil((4E-05*fck**2 - 0.0071*fck + 0.634)*100)/100
    return [min(alpha,0.75), min(beta, 0.38)]

def EC2_lambda_eta(fck, units="MPa"):
    """ [lambda, eta] = EC2_lambda_eta(fck) """
    fck = convert_2_MPa(fck, units)
    if fck > 90:
        raise ValueError('fck = {} {} > 90 Mpa'.format(fck, units))
    elif fck > 50:
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
#   Miscellaneous
#===========================================================================


def convert_2_psi(fck, units):
    if units == "MPa":
        return 145.038*fck
    elif units == "psi":
        return fck
    else:
        raise KeyError("Only psi or MPa")

def convert_2_pcf(rho, units):
    if units == "psi":
        return rho
    elif units == "MPa":
        return 6.366*rho
    else:
        raise KeyError("Only psi or MPa")

def convert_2_MPa(fck, units):
    if units == "MPa":
        return fck
    elif units == "psi":
        return fck/145.038
    else:
        raise KeyError("Only psi or MPa")

def convert_2_mm(h, units):
    if units == "MPa":
        return h
    elif units == "psi":
        return 25.4*h
    else:
        raise KeyError("Only psi or MPa")

if __name__ == "__main__":
    print("\n Ec for C30/37: \n")
    print(ACI_elastic_modulus(30, units="MPa"))
    print(EC2_elastic_modulus(30))
    print("\n Tensile strength for C30/37: \n")
    print(ACI_tensile_strength(30,units="MPa"))
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

