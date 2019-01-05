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
		gamma_S = 1.15
		gamma_P = 1.1
		niter = 20

	Instance Members:
	
	Public Methods:
		[alpha, beta, la, eta] = design_parameters(self)
		[fcd, fyd, fp01d, fpd] = design_properties(self)
	
	Static Methods:
		plotArrow(self, ax, x, y, dx, dy, color)
		plotStress(self, ax, fcd, h_top, h_bot) 
	"""
	gamma_P = 1.1
	Pm = Ap = 0

	def __init__(self):
		ReinforcedBeam.__init__(self)
		self.fig = plt.figure(figsize=(9, 3))
		self.sigma_p = self.Pm/self.Ap
		self.f = self.d_P - self.d_P0

	def design_parameters(self):
		[alpha, beta] = EC2.alpha_beta(self.fck, self.units)
		[la, eta] = EC2.lambda_eta(self.fck, self.units)
		logging.debug("    a = {:3.2f}, b = {:3.2f}".format(alpha, beta))
		logging.debug("    lambda = {:3.2f}, eta = {:3.2f}".format(la, eta))
		return [alpha, beta, la, eta]

	def design_properties(self):
		fcd = self.fck/self.gamma_c
		fyd = self.fyk/self.gamma_S
		fp01d = 0.9*self.fpk/self.gamma_P
		fpd = self.fpk/self.gamma_P
		logging.debug("    fcd = {:6.2f}, fyd = {:6.2f}".format(fcd, fyd))
		logging.debug("    fp01d = {:6.2f}, fpd = {:6.2f}".format(fp01d, fpd))
		return [fcd, fyd, fp01d, fpd]

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
		[alpha, beta, la, eta] = self.design_parameters()
		[fcd, fyd, fp01d, fpd] = self.design_properties()
		ecu3 = EC2.ultimate_strain(self.fck, self.units)
		epy = fp01d/self.Ep
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
