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
    [fck, fyd] = convert_2_MPa(np.array([fck, fyd]), units)
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
