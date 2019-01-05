# Import Standard Libraries
import logging
import scipy as np
import matplotlib.pyplot as plt

# Import Local Libraries
import Util_ACI as ACI
import Util_EC2 as EC2
from PSC import PrestressedBeam
from TeeBeam import TeeBeam

#===========================================================================
#   BoxGirder
#===========================================================================


class BoxGirder(PrestressedBeam):
	Bt = 1500
	ht = 250

	Bb = 1500
	hb = 250
	
	bt = 1000
	bb = 1000
	H = 1000
	
	d_S = 800
	d_P = 800
	d_P0 = 300
	As = 0
	Ap = 10080

	fck = 50
	fyk = 500
	fpk = 1860
	
	Es = 205000 
	k = 0.01
	Ep = 195000
	Pm = 11344030
	units = "MPa"

	#   |<--------------- Bt --------------->|
	#	______________________________________  ___
	#	|                       |            |   |
	#	|                       |            |   |   ___
	#	|                       ht           |   |    |
	#	|         ______________|___         |   |    |
	#   |<-bwt/2->|<------bt----|->|<-bwt/2->|   |    |
	#	|         |             |  |         |   |    |
	#   |         |            hw  |         |   H    |   ___
	#   |         |             |  |         |   |   yft   |
	#   |         |<------bb----|->|         |   |    |    |
	#	|<-bwb/2->|_____________|__|<-bwb/2->|   |    |   yfw
	#	|                       |            |   |    |    |   ___
	#	|                       |            |   |    |    |    |
	#	|                       hb           |   |    |    |   yfb
	#	|_______________________|____________|  _|_  _|_  _|_  _|_
	#
	#   |<--------------- Bb --------------->|

	def __init__(self):

		PrestressedBeam.__init__(self)

		# Slope of girder width
		slope = (self.Bb - self.Bt)/self.H		
		
		# Properties of top flange
		Bt2 = self.Bt + slope*self.ht
		Aft = (self.Bt + Bt2)*self.ht/2
		Ift = self.momentInertia(self.ht, self.Bt, Bt2)
		yft = self.H - self.neutralAxis(self.ht, self.Bt, Bt2)
		logging.debug("    Aft = {:6.2f}, Ift = {:10.2f}".format(Aft, Ift))

		# Properties of bottom flange
		Bb2 = self.Bb - slope*self.hb
		Afb = (self.Bb + Bb2)*self.hb/2
		Ifb = self.momentInertia(self.hb, self.Bb, Bb2)
		yfb = self.neutralAxis(self.hb, self.Bb, Bb2)
		logging.debug("    Afb = {:6.2f}, Ifb = {:10.2f}".format(Afb, Ifb))

		# Properties of web
		hw = self.H - self.ht - self.hb
		bwt = Bt2 - self.bt
		bwb = Bb2 - self.bb
		Aw = (bwt + bwb)*hw/2
		Iw = self.momentInertia(hw, bwb, bwt)
		yw = self.hb + self.neutralAxis(hw, bwb, bwt)
		logging.debug("    Aw = {:6.2f}, Iw = {:10.2f}".format(Aw, Iw))

		# Properties of BoxGirder
		self.Ac = Aft + Aw + Afb
		self.yb = (Aft*yft + Aw*yw + Afb*yfb)/self.Ac
		self.yt = self.H - self.yb
		logging.debug("    yt = {:6.2f}, yb = {:}".format(self.yb, self.yt))

		self.e_P = self.d_P - self.yt
		logging.debug("    e_P = {:6.2f}".format(self.e_P))

		self.I = Ift + Iw + Ifb + Aft*(yft - self.yb)**2 \
			+ Aw*(yw - self.yb)**2 + Afb*(yfb - self.yb)**2
		logging.debug("    Ac = {:6.2f}, I = {:10.2f}".format(self.Ac, self.I))

		self.Wt = self.I/self.yt
		self.Wb = self.I/self.yb

		logging.debug("    Wt = {:6.2f}, Wb = {:6.2f}".format(self.Wt, self.Wb))

	def EC2_design_moment(self):
		[alpha, beta, la, eta] = self.design_parameters()
		[fcd, fyd, fp01d, fpd] = self.design_properties()
		ecu3 = EC2.ultimate_strain(self.fck, self.units)
		epy = fp01d/self.Ep
		epu = 0.035
		ey = fyd/self.Es
		ax = self.fig.add_subplot(132)
		ax.plot([0, 0], [0, self.h], linewidth=2, color='k')
		# Assume fp = 0.95*fpd and fs = fyd
		fp = 0.95*fpd
		fs = fyd
		logging.debug("    Assume fp = 0.95 fpd = {:6.2f}".format(fp))

		Xu = (self.Ap*(fp-self.sigma_p)+self.Pm)/(alpha*self.bf*fcd)

		for iter in range(self.niter):

			if Xu < self.ht:
				Xu = (self.Ap*(fp-self.sigma_p)+self.Pm)/(alpha*self.bf*fcd)

			elif Xu > self.ht:
				X = (self.Ap*fp-eta*fcd*self.ht *
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
		if Xu < self.ht:
			# _________
			#   |  |     ||
			#   |  |     ||   b*Xu
			#   |  |_    ||<------ Nc
			#   |    \_  ||    |
			#   |      \_||    y1
			#   |        ||    |
			#   dp ----- ||<------ Pm　＠　neutral axis
			#   |    |   ||    y2
			#   |    e   ||<------ Ns
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

		elif Xu > self.ht:
			# _________
			#   |  |     ||
			#   |  |     ||<----------- Nc1
			#   |  |_____||<---- Nc2 |
			#   |  |     ||   |      |
			#   |  |_____||   y2     y1
			#   |        ||   |      |
			#   dp ----  ||<------------ Pm　＠　neutral axis
			#   |    |   ||    y3
			#   |    e   ||<------ Ns
			#   |    |   ||
			#  --------- ||<------Delta Pm
			#            ||

			Nc1 = (self.bf-self.bw)*self.ht*fcd
			y1 = self.d_P - self.e_P - self.ht/2
			self.plotStress(ax, fcd, self.h, self.h - self.ht)
			self.plotArrow(ax, Nc1/10**5, self.yb+y1, -Nc1/10**5, 0, 'b')
			logging.info("    Nc1 = {:6.2f}, y1 = {:6.2f}".format(Nc1, y1))

			Nc2 = self.bw*X*fcd
			y2 = self.d_P - self.e_P - X/2
			self.plotStress(ax, fcd, self.h, self.h - Xu)
			self.plotArrow(ax, Nc2/10**5, self.yb+y2, -Nc2/10**5, 0, 'b')
			logging.info("    Nc2 = {:6.2f}, y2 = {:6.2f}".format(Nc2, y2))

			Ns = self.As*fs
			y3 = self.d_P - y1 - self.ht/2
			logging.info("    Ns = {:6.2f}, y3 = {:6.2f}".format(Ns, y3))

			Delta_P = self.Ap*fp - self.Pm
			self.plotArrow(ax, 0, self.h - self.d_P, Delta_P/10**5, 0, 'r')
			logging.info(
				"    Delta_P = {:6.2f}, e_P = {:6.2f}".format(Delta_P, self.e_P))

			MRd = Nc1*y1 + Nc2*y2 + Ns*y3 + Delta_P*self.e_P
			logging.info("    MRd = {:5.2f}".format(MRd))
			plt.xlim((-1.2*fcd, 2*fcd))
		return MRd

	def plotBoxGirder(self):
		x0 = self.bf
		x1 = (self.bf + self.bw)/2
		x2 = (self.bf - self.bw)/2

		y0 = self.h
		y1 = self.h - self.ht
		X = [0, x0, x0, x1, x1, x2, x2, 0, 0]
		Y = [y0, y0, y1, y1, 0, 0, y1, y1, y0]

		ax = self.fig.add_subplot(131)
		ax.plot(X, Y, linewidth=1, color='k')
		ax.fill_between(X, Y, alpha=0.5)

		ax.plot([0, self.bf], [self.yb, self.yb], linestyle='--', color='k')

		ax.plot(self.bf/2, self.yb, marker='+', color='k')

		if self.Ap != 0:
			ax.plot(self.bf/2, self.h-self.d_P, marker='o', color='r')

		if self.As != 0:
			ax.plot(self.bf/2, self.h-self.d_S, marker='o', color='k')

		ax.plot(self.bf/2, self.h-self.d_P0, marker='o', color='r')
		plt.show()

	def plotStressDist(self, Xu, *args):

		ax = self.fig.add_subplot(132)
		ax.plot([0, 0], [0, self.h], linewidth=2, color='k')

		# Moment Capacity
		if Xu < self.ht:
			Nc1 = args[0]
			print(Nc1)
		elif Xu > self.ht:

			[fcd, Nc1, y1, Nc2, y2, Ns, y3] = args

			h0 = self.h - self.ht
			self.__plotConcreteStress(ax, fcd, self.h, h0)

			h1 = self.h - Xu
			self.__plotConcreteStress(ax, fcd, self.h, h1)

			# self.__plotSteelStress(ax, Ns, y3)
			print(Ns)
			# y1 = [self.h, self.h, h1, h1, self.h]

			h1 = h1 + Xu/2

			plt.xlim((-1.2*fcd, 1.2*fcd))
			plt.show()


if __name__ == "__main__":
	logging.basicConfig(level=logging.DEBUG, format='%(message)s')
	beam = BoxGirder()
	beam.EC2_design_moment()
	beam.plotBoxGirder()
	# beam.plotStressDist()
