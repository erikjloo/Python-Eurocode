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

		logging.debug("    Wt = {:6.2f}, Wb = {:6.2f}".format(self.Wt, self.Wb))

	def plotTeeBeam(self):
		x1 = self.bf
		x2 = (self.bf + self.bw)/2
		x3 = (self.bf - self.bw)/2

		y0 = self.H
		y1 = self.H - self.hf
		X = [0, x1, x1, x2, x2, x3, x3, 0, 0]
		Y = [y0, y0, y1, y1, 0, 0, y1, y1, y0]

		ax = self.fig.add_subplot(131)
		ax.plot(X, Y, linewidth=1, color='k')
		ax.fill(X, Y, color = 'lightblue')

		ax.plot([0, self.bf], [self.yb, self.yb], linestyle='--', color='k')

		ax.plot(self.bf/2, self.yb, marker='+', color='k')

		if self.Ap != 0:
			ax.plot(self.bf/2, self.H - self.d_P, marker='o', color='r')
			ax.plot(self.bf/2, self.H - self.d_P0, marker='o', color='r')

		if self.As != 0:
			ax.plot(self.bf/2, self.H - self.d_S, marker='o', color='k')

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
