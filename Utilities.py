#===========================================================================
#   Utilities
#===========================================================================


def convert_2_psi(fck, units):
    if units not in ["MPa", "psi"]:
        raise KeyError("Only psi or MPa")
    return fck if units == "psi" else 145.038*fck


def convert_2_in(h, units):
    if units not in ["MPa", "psi"]:
        raise KeyError("Only psi or MPa")
    return h if units == "psi" else h/25.4


def convert_2_in2(A, units):
    if units not in ["MPa", "psi"]:
        raise KeyError("Only psi or MPa")
    return A if units == "psi" else A/(25.4**2)


def convert_2_pcf(rho, units):
    if units not in ["MPa", "psi"]:
        raise KeyError("Only psi or MPa")
    return rho if units == "psi" else 6.366*rho


def convert_2_MPa(fck, units):
    if units not in ["MPa", "psi"]:
        raise KeyError("Only psi or MPa")
    return fck if units == "MPa" else fck/145.038


def convert_2_mm(h, units):
    if units not in ["MPa", "psi"]:
        raise KeyError("Only psi or MPa")
    return h if units == "MPa" else 25.4*h


def convert_2_mm2(A, units):
    if units not in ["MPa", "psi"]:
        raise KeyError("Only psi or MPa")
    return A if units == "MPa" else 25.4**2*A
