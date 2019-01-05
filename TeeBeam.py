# Import Standard Libraries
import logging
import scipy as np
import matplotlib.pyplot as plt

# Import Local Libraries
import Util_ACI as ACI
import Util_EC2 as EC2
from PSC import PrestressedBeam

#===========================================================================
#   T Beam
#===========================================================================


class TeeBeam(PrestressedBeam):
	bf = 1200
	hf = 150

	bw = 300
	H = 1200

	d_S = 1000
	d_P = 1000
	d_P0 = 500

	As = 0
	Ap = 2400

	fck = 30
	fyk = 500
	fpk = 1860

	Es = 205000
	k = 0.01
	Ep = 195000
	Pm = 2772000
	units = "MPa"

    #   |<------- bf ------->|
    #   ______________________  ___
    #   |                |   |   |
    #   |               hf   |   |   ___
    #   |_______      ___|___|   |    |
    #          |      |  |       |    |
    #          |<-bw->|  |       |    |
    #          |      |  |       H    |   ___
    #          |      |  hw      |    yf   |
    #          |      |  |       |    |    yw
    #          |      |  |       |    |    |
    #          |______| _|_     _|_  _|_  _|_

	def __init__(self):

		PrestressedBeam.__init__(self)

        # Properties of flange
		Af = self.hf*self.bf
		If = self.bf*self.hf**3/12
		yf = self.H - self.hf/2
		logging.debug("    Af = {:6.2f}, If = {:10.2f}".format(Af,If))

		# Properties of web
		hw = self.H - self.hf
		Aw = self.bw*hw
		Iw = self.bw*hw**3/12
		yw = hw/2
		logging.debug("    Aw = {:6.2f}, Iw = {:10.2f}".format(Aw, Iw))

		# Properties of TeeBeam
		self.Ac = Af + Aw
		self.yb = (Af*yf+Aw*yw)/self.Ac
		self.yt = self.H - self.yb
		logging.debug("    yt = {:6.2f}, yb = {:}".format(self.yb, self.yt))

		self.e_P = self.d_P - self.yt
		logging.debug("    e_P = {:6.2f}".format(self.e_P))

		self.I = Iw + If + Af*(self.yb - yf)**2 + Aw*(self.yb - yw)**2
		logging.debug("    Ac = {:6.2f}, I = {:10.2f}".format(self.Ac, self.I))

		self.Wt = self.I/self.yt
		self.Wb = self.I/self.yb

		logging.debug(
			"    Wt = {:6.2f}, Wb = {:6.2f}".format(self.Wt, self.Wb))

	def EC2_design_moment(self):
		[alpha, beta, la, eta] = self.design_parameters()
		[fcd, fyd, fp01d, fpd] = self.design_properties()
		ecu3 = EC2.ultimate_strain(self.fck, self.units)
		epy = fp01d/self.Ep
		epu = 0.035
		ey = fyd/self.Es
		ax = self.fig.add_subplot(132)
		ax.plot([0, 0], [0, self.H], linewidth=2, color='k')
		# Assume fp = 0.95*fpd and fs = fyd
		fp = 0.95*fpd
		fs = fyd
		logging.debug("    Assume fp = 0.95 fpd = {:6.2f}".format(fp))

		Xu = (self.Ap*(fp-self.sigma_p)+self.Pm)/(alpha*self.bf*fcd)

		for iter in range(self.niter):

			if Xu < self.hf:
				Xu = (self.Ap*(fp-self.sigma_p)+self.Pm)/(alpha*self.bf*fcd)

			elif Xu > self.hf:
				X = (self.Ap*fp-eta*fcd*self.hf *
                                    (self.bf-self.bw))/(eta*fcd*self.bw)
				Xu = X/la

			logging.debug(
				"Iteration {:1.0f}: \n    Xu = {:6.2f}".format(iter, Xu))

			# Is the assumption correct?
			Delta_ep = (self.d_P-Xu)*ecu3/Xu
			ep = self.sigma_p/self.Ep + Delta_ep
			fp_new = fp01d + (ep-epy)*(fpd-fp01d)/(epu-epy)
			logging.debug("    ep = {:6.5f}, fp = {:6.2f}".format(ep, fp_new))

			es = (self.d_S-Xu)*ecu3/Xu
			fs_new = fyd + (es-ey)*self.k

			# Approximation Error

			if abs((fp_new-fp)/fp_new) < 0.001 and abs((fs_new-fs)/fs_new) < 0.001:
				break

			# New assumption f = f_new
			fp = fp_new
			fs = fs_new

		# Moment Capacity
		if Xu < self.hf:
			# _________
			#   |  |     ||
			#   |  |     ||   b*Xu
			#   |  |_    ||<------ Nc
			#   |    \_  ||    |
			#   |      \_||    y1
			#   |        ||    |
			#   dp ----- ||<------ Pm　＠　neutral axis
			#   |    |   ||    y2
			#   |   e_P  ||<------ Ns
			#   |    |   ||
			#  --------- ||<------Delta Pm
			#            ||

			Nc = alpha*self.bf*Xu*fcd
			y1 = self.d_P - self.e_P - beta*Xu
			logging.info("    Nc = {:6.2f}, y1 = {:6.2f}".format(Nc, y1))

			Ns = self.As*fs
			y2 = self.d_S - y1 - beta*Xu
			logging.info("    Ns = {:6.2f}, y2 = {:6.2f}".format(Ns, y2))

			Delta_P = self.Ap*fp - self.Pm
			logging.info(
				"    Delta_P = {:6.2f}, e_P = {:6.2f}".format(Delta_P, self.e_P))

			MRd = Nc*y1 + Ns*y2 + Delta_P*self.e_P

		elif Xu > self.hf:
			# _________
			#   |  |     ||
			#   |  |     ||<----------- Nc1
			#   |  |_____||<---- Nc2 |
			#   |  |     ||   |      |
			#   |  |_____||   y2     y1
			#   |        ||   |      |
			#   dp ----  ||<------------ Pm　＠　neutral axis
			#   |    |   ||    y3
			#   |   e_P  ||<------ Ns
			#   |    |   ||
			#  --------- ||<------Delta Pm
			#            ||

			Nc1 = (self.bf-self.bw)*self.hf*fcd
			y1 = self.d_P - self.e_P - self.hf/2
			self.plotStress(ax, fcd, self.H, self.H - self.hf)
			self.plotArrow(ax, Nc1/10**5, self.yb+y1, -Nc1/10**5, 0, 'b')
			logging.info("    Nc1 = {:6.2f}, y1 = {:6.2f}".format(Nc1, y1))

			Nc2 = self.bw*X*fcd
			y2 = self.d_P - self.e_P - X/2
			self.plotStress(ax, fcd, self.H, self.H - Xu)
			self.plotArrow(ax, Nc2/10**5, self.yb+y2, -Nc2/10**5, 0, 'b')
			logging.info("    Nc2 = {:6.2f}, y2 = {:6.2f}".format(Nc2, y2))

			Ns = self.As*fs
			y3 = self.d_P - y1 - self.hf/2
			logging.info("    Ns = {:6.2f}, y3 = {:6.2f}".format(Ns, y3))

			Delta_P = self.Ap*fp - self.Pm
			self.plotArrow(ax, 0, self.H - self.d_P, Delta_P/10**5, 0, 'r')
			logging.info(
				"    Delta_P = {:6.2f}, e_P = {:6.2f}".format(Delta_P, self.e_P))

			MRd = Nc1*y1 + Nc2*y2 + Ns*y3 + Delta_P*self.e_P
			logging.info("    MRd = {:5.2f}".format(MRd))
			plt.xlim((-1.2*fcd, 2*fcd))
		return MRd

	def plotTeeBeam(self):
		x0 = self.bf
		x1 = (self.bf + self.bw)/2
		x2 = (self.bf - self.bw)/2

		y0 = self.H
		y1 = self.H - self.hf
		X = [0, x0, x0, x1, x1, x2, x2, 0, 0]
		Y = [y0, y0, y1, y1, 0, 0, y1, y1, y0]

		ax = self.fig.add_subplot(131)
		ax.plot(X, Y, linewidth=1, color='k')
		ax.fill_between(X, Y, alpha=0.5)

		ax.plot([0, self.bf], [self.yb, self.yb], linestyle='--', color='k')

		ax.plot(self.bf/2, self.yb, marker='+', color='k')

		if self.Ap != 0:
			ax.plot(self.bf/2, self.H - self.d_P, marker='o', color='r')

		if self.As != 0:
			ax.plot(self.bf/2, self.H - self.d_S, marker='o', color='k')

		ax.plot(self.bf/2, self.H - self.d_P0, marker='o', color='r')
		plt.show()

	def plotStressDist(self, Xu, *args):

		ax = self.fig.add_subplot(132)
		ax.plot([0, 0], [0, self.H], linewidth=2, color='k')

		# Moment Capacity
		if Xu < self.hf:
			Nc1 = args[0]
			print(Nc1)
		elif Xu > self.hf:

			[fcd, Nc1, y1, Nc2, y2, Ns, y3] = args

			h0 = self.H - self.hf
			self.__plotConcreteStress(ax, fcd, self.H, h0)

			h1 = self.H - Xu
			self.__plotConcreteStress(ax, fcd, self.H, h1)

			print(Ns)

			h1 = h1 + Xu/2

			plt.xlim((-1.2*fcd, 1.2*fcd))
			plt.show()


if __name__ == "__main__":
	logging.basicConfig(level=logging.DEBUG, format='%(message)s')
	beam = TeeBeam()
	beam.EC2_design_moment()
	beam.plotTeeBeam()
	# beam.plotStressDist()
