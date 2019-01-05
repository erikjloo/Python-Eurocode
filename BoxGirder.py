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

	Bb = 1200
	hb = 250
	
	bt = 1000
	bb = 800
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
		self.slope = (self.Bt - self.Bb)/self.H		
		
		# Properties of top flange
		Bt2 = self.Bt - self.slope*self.ht
		self.hf = self.ht # Needed for EC2_design_moment
		self.bf = (self.Bt + Bt2)/2
		Aft = (self.Bt + Bt2)*self.ht/2
		Ift = self.momentInertia(self.ht, self.Bt, Bt2)
		yft = self.H - self.neutralAxis(self.ht, self.Bt, Bt2)
		logging.debug("    Aft = {:6.2f}, Ift = {:10.2f}".format(Aft, Ift))

		# Properties of bottom flange
		Bb2 = self.Bb + self.slope*self.hb
		Afb = (self.Bb + Bb2)*self.hb/2
		Ifb = self.momentInertia(self.hb, self.Bb, Bb2)
		yfb = self.neutralAxis(self.hb, self.Bb, Bb2)
		logging.debug("    Afb = {:6.2f}, Ifb = {:10.2f}".format(Afb, Ifb))

		# Properties of web
		hw = self.H - self.ht - self.hb
		self.bwt = Bt2 - self.bt
		self.bwb = Bb2 - self.bb
		self.bw = (self.bwt + self.bwb)/2 # Needed for EC2_design_moment
		Aw = (self.bwt + self.bwb)*hw/2
		Iw = self.momentInertia(hw, self.bwb, self.bwt)
		yw = self.hb + self.neutralAxis(hw, self.bwb, self.bwt)
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

	def plotBoxGirder(self):
		X1 = self.Bt
		X3 = self.H*self.slope/2
		X2 = X3 + self.Bb

		Y0 = self.H

		x0 = self.ht*self.slope/2 + self.bwt/2
		x1 = x0 + self.bt
		x3 = (self.H - self.hb)*self.slope/2 + self.bwb/2
		x2 = x3 + self.bb

		y0 = self.H - self.ht
		y1 = self.hb

		ax = self.fig.add_subplot(131)
		X = [0, X1, X2, X3, 0]
		Y = [Y0, Y0, 0, 0, Y0]
		ax.plot(X, Y, linewidth=1, color='k')
		ax.fill(X, Y, color = 'lightblue')

		X = [x0, x1, x2, x3, x0]
		Y = [y0, y0, y1, y1, y0]
		ax.plot(X, Y, linewidth=1, color='k')
		ax.fill(X, Y, color = 'w')

		ax.plot([0, self.Bt], [self.yb, self.yb], linestyle='--', color='k')

		ax.plot(self.Bt/2, self.yb, marker='+', color='k')

		if self.Ap != 0:
			ax.plot(self.Bt/2, self.H-self.d_P, marker='o', color='r')
			ax.plot(self.Bt/2, self.H-self.d_P0, marker='o', color='r')

		if self.As != 0:
			ax.plot(self.Bt/2, self.H-self.d_S, marker='o', color='k')

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
