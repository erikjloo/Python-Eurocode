# Import Standard Libraries
from abc import ABCMeta, abstractmethod
import logging
import scipy as np

# Import Local Libraries
import Util_ACI as ACI
import Util_EC2 as EC2


#===========================================================================
#   Reinforced Beam
#===========================================================================

class ReinforcedBeam(metaclass=ABCMeta):
	""" Abstract Reinforced Beam Class
	
	Static Members:
		gamma_c = 1.5
		gamma_S = 1.15
		niter = 20

	Instance Members:
	
	Static Methods:
		plotArrow(self, ax, x, y, dx, dy, color)
		plotStress(self, ax, fcd, h_top, h_bot)
	
	"""

	gamma_c = 1.5
	gamma_S = 1.15
	niter = 20

	H = d_S = As = fck = fyk = Es = k = 0
	units = "MPa"

	def __init__(self, name=None):
		self.ta = self.H/500
		self.name = name

	def plotArrow(self, ax, x, y, dx, dy, color):
		if x == 0:
			ax.arrow(x, y, dx - self.ta, dy, fc=color, ec=color,
                            width=10*self.ta, head_width=10*self.ta, head_length=self.ta)
		else:
			ax.arrow(x, y, dx + self.ta, dy, fc=color, ec=color,
                            width=self.ta, head_width=10*self.ta, head_length=self.ta)

	@staticmethod
	def plotStress(ax, fcd, h_top, h_bot):
		X0 = [-fcd, 0]
		Y0 = [h_top, h_top]
		Y1 = [h_bot, h_bot]
		ax.fill_between(X0, Y0, Y1, color='b', alpha=0.3)

	@staticmethod
	def neutralAxis(h, bb, bt):
		""" Input:	h = trapezoidal height
					bb = bottom width
					bt = top width
			Output:	yc = bottom to neutral axis """
		return h*(2*bt + bb)/3/(bb + bt)

	@staticmethod
	def momentInertia(h, bb, bt):
		""" Input:	h = trapezoidal height
					bb = bottom width
					bt = top width
			Output:	I = moment of inertia about neutral axis """
		return h**3*(bt**2+4*bt*bb+bb**2)/36/(bt+bb)

	@staticmethod
	def beamFactory():
		pass

#===========================================================================
#   Rectangular Beam
#===========================================================================


class RectangularBeam(ReinforcedBeam):
	# b = 558.8
	# h = 609.6
	# d = 546.1
	# As = 3870
	# fck = 27.6
	# fyk = 414
	# rho = 25
	# Es = 200000
	# units = "MPa"
	b = 22
	h = 24
	d = 21.5
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

		logging.debug("Uncracked moment capacity per ACI")
		Ec = ACI.elastic_modulus(self.fck, self.rho, self.units)
		fr = ACI.tensile_strength(self.fck, self.units)
		n = self.Es/Ec
		logging.debug("    fr = {:6.2f}, Es/Ec = {:3.2f}".format(fr, n))

		Ac = self.b*self.h
		As = (n-1)*self.As
		logging.debug("    b*h = {:6.2f}, (n-1)*As = {:6.2f}".format(Ac, As))

		y_bott = (Ac*self.h/2 + As*(self.h-self.d))/(Ac+As)
		I_uncr = self.b*self.h**3/12\
			+ Ac*(self.h/2-y_bott)**2\
			+ As*(y_bott-(self.h-self.d))**2
		logging.debug("    Botton to NA: y_bott = {:6.2f}".format(y_bott))
		logging.debug("    Moment of inertia = {:10.2f} ".format(I_uncr))

		Mcr = fr*I_uncr/y_bott
		logging.info("    Mcr = {:5.2f}".format(Mcr))
		return Mcr

	def ACI_elastic_moment(self):

		Ec = ACI.elastic_modulus(self.fck, self.rho, self.units)
		fc = 0.5*self.fck
		n = self.Es/Ec
		As = n*self.As
		fs = self.fyk/n
		logging.debug("    n = {:3.2f}, n*As = {:6.2f}".format(n, As))

		kd = (-As+np.sqrt(As**2+2*self.b*As*self.d))/(self.b)
		Icr = self.b*kd**3/12 + self.b*kd*(kd/2)**2 + As*(self.d-kd)**2
		logging.debug("    Top to NA: kd = {:6.2f}".format(kd))
		logging.debug("    Moment of inertia = {:10.2f}".format(Icr))

		Mel = min(fc*Icr/kd, fs*Icr/(self.d-kd))
		logging.info("    Mel = {:5.2f}".format(Mel))
		return Mel

	def ACI_design_moment(self):

		beta1 = ACI.beta(self.fck, self.units)
		logging.debug("    beta1 = {:3.2f}".format(beta1))
		c = (self.As*self.fyk)/(0.85*self.fck*self.b)/beta1
		logging.debug("    c = {:6.2f}".format(c))
		phi = ACI.ductility_requirement(c, self.d, type="beam")
		logging.debug("    phi = {:3.2f}".format(phi))
		ACI.steel_ratio(self.As, self.fck, self.fyk,
                  self.b, self.d, self.units)
		MRd = phi*self.As*self.fyk*(self.d-beta1*c/2)
		logging.info("    MRd = {:5.2f}".format(MRd))
		return MRd

	#---------------------------------------------------------------------------
	#   EC2 Equations
	#---------------------------------------------------------------------------

	def EC2_cracking_moment(self):

		logging.debug("Uncracked moment capacity per EC2")
		Ec = EC2.elastic_modulus(self.fck, self.units)
		fr = EC2.flex_tensile_strength(self.fck, self.h, self.units)
		n = self.Es/Ec
		logging.debug("    fr = {:6.2f}, Es/Ec = {:3.2f}".format(fr, n))

		Ac = self.b*self.h
		As = (n-1)*self.As
		logging.debug("    b*h = {:6.2f}, (n-1)*As = {:6.2f}".format(Ac, As))

		y_bott = (Ac*self.h/2 + As*(self.h-self.d))/(Ac+As)
		I_uncr = self.b*self.h**3/12\
			+ Ac*(self.h/2-y_bott)**2\
			+ As*(y_bott-(self.h-self.d))**2
		logging.debug("    Botton to NA: y_bott = {:6.2f}".format(y_bott))
		logging.debug("    Moment of inertia = {:10.2f} ".format(I_uncr))

		Mcr = fr*I_uncr/y_bott
		logging.info("    Mcr = {:5.2f}".format(Mcr))
		return Mcr

	def EC2_elastic_moment(self):

		Ec = EC2.elastic_modulus(self.fck, self.units)
		fc = 0.5*self.fck
		n = self.Es/Ec
		As = n*self.As
		fs = self.fyk/n
		logging.debug("    n = {:3.2f}, n*As = {:6.2f}".format(n, As))

		kd = (-As+np.sqrt(As**2+2*self.b*As*self.d))/(self.b)
		Icr = self.b*kd**3/12 + self.b*kd*(kd/2)**2 + As*(self.d-kd)**2
		logging.debug("    Top to NA: kd = {:6.2f}".format(kd))
		logging.debug("    Moment of inertia = {:10.2f}".format(Icr))

		Mel = min(fc*Icr/kd, fs*Icr/(self.d-kd))
		logging.info("    Mel = {:5.2f}".format(Mel))
		return Mel

	def EC2_design_moment(self):
		[alpha, beta] = EC2.alpha_beta(self.fck, self.units)
		logging.debug("    a = {:3.2f}, b = {:3.2f}".format(alpha, beta))
		fcd = self.fck/self.gamma_c
		fyd = self.fyk/self.gamma_S
		logging.debug("    fcd = {:6.2f}, fyd = {:6.2f}".format(fcd, fyd))
		Xu = self.As*fyd/(alpha*self.b*fcd)
		Xu_max = EC2.ductility_requirement(
			Xu, self.d, self.fck, fyd, self.units)

		EC2.steel_ratio(self.As, self.fck, self.fyk, self.b,
                  self.d, self.h, Xu_max, self.units)

		MRd = self.As*fyd*(self.d-beta*Xu)
		logging.info("    MRd = {:5.2f}".format(MRd))
		return MRd


class DoublyReinforcedBeam():
	a = 2
	h = 609.6
	d_1 = 20
	As_1 = 10
	d_2 = 546.1
	As_2 = 40
	b = 558.8
	As = 2.8575**2*np.pi/4
	fck = 27.6
	fyk = 414
	rho = 25
	Es = 200000
	units = "MPa"


if __name__ == "__main__":
	logging.basicConfig(level=logging.DEBUG, format='%(message)s')

	# print("\n Ec for C30/37: \n")
	# print(ACI.elastic_modulus(30, rho=25, units="MPa"))
	# print(EC2.elastic_modulus(30))
	# print("\n Tensile strength for C30/37: \n")
	# print(ACI.tensile_strength(30, units="MPa"))
	# print(EC2.tensile_strength(30))
	# print("\n Tensile strength for C55/67: \n")
	# print(ACI.tensile_strength(55, units="MPa"))
	# print(EC2.tensile_strength(55))
	# print("\n Ultimate strain for C30/37 \n")
	# print(ACI.ultimate_strain(30, units="MPa"))
	# print(EC2.ultimate_strain(30))
	# print("\n Ultimate strain for C55/67 \n")
	# print(ACI.ultimate_strain(55, units="MPa"))
	# print(EC2.ultimate_strain(55))
	# print("\n alpha & Beta factors: \n")
	# print(EC2.alpha_beta(10))
	# print(EC2.alpha_beta(60))
	# print(EC2.alpha_beta(80))
	# print("\n Lambda & Eta factors: \n")
	# print(EC2.lambda_eta(10))
	# print(EC2.lambda_eta(60))
	# print(EC2.lambda_eta(80))

	beam = RectangularBeam()
	# print('\n Cracking Moment: \n')
	# MRd = beam.ACI.cracking_moment()
	# print(MRd/12000)
	# MRd = beam.EC2.cracking_moment()
	# print(MRd/12000)

	print('\n Elastic Moment: \n')
	MRd = beam.ACI_elastic_moment()
	print(MRd/12000)
	MRd = beam.EC2_elastic_moment()
	print(MRd/12000)

	print('\n Design Moment: \n')
	MRd = beam.ACI_design_moment()
	print(MRd/12000)
	MRd = beam.EC2_design_moment()
	print(MRd/12000)
