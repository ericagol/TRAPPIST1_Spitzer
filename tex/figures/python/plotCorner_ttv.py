import matplotlib
#matplotlib.use('PS')
#matplotlib.use('agg')

import pylab as pl
import numpy as np
import math

pl.rc('font', size=30)
params = {'legend.fontsize': 30}
pl.rcParams.update(params)
pl.rc('text', usetex=True)

#yt = np.arange(1.0)*0.1
#xt = np.arange(4.0, 7.1, 1.0)
pl.figure(figsize=(30, 30))


#file = 'MCMCaS.dat'
#file = 'state_total.txt'
file = '../../../data/state_total.txt'

n = 35

for i in range(0, n, 1):
	for j in range(0, i + 1, 1):
	#for j in range(i, i + 1, 1):

#		pi, pj, chi =  np.loadtxt(file, unpack=True, usecols = (i, j, 35))
		pi, pj =  np.loadtxt(file, unpack=True, usecols = (i, j))


		t = []
		t.append("b")
		t.append("c")
		t.append("d")
		t.append("e")
		t.append("f")
		t.append("g")
		t.append("h")


		if(j % 5 == 0):
			xl = r"$m_%s$" % t[j // 5]
		if(j % 5 == 1):
			xl = r"$P_%s$" % t[j // 5]
		if(j % 5 == 2):
			xl = r"$t0_%s$" % t[j // 5]
		if(j % 5 == 3):
			xl = r"$k_%s$" % t[j // 5]
		if(j % 5 == 4):
			xl = r"$h_%s$" % t[j // 5]


		if(i % 5 == 0):
			yl = r"$m_%s$" % t[i // 5]
		if(i % 5 == 1):
			yl = r"$P_%s$" % t[i // 5]
		if(i % 5 == 2):
			yl = r"$t0_%s$" % t[i // 5]
		if(i % 5 == 3):
			yl = r"$k_%s$" % t[i // 5]
		if(i % 5 == 4):
			yl = r"$h_%s$" % t[i // 5]
		
		

		i1 = np.max(pi)
		i0 = np.min(pi)
		j1 = np.max(pj)
		j0 = np.min(pj)

		print(i, j, i * n + j)

		ax1=pl.subplot(n,n,i * n + j + 1)
		col = 'b'

		if(i != j):
						
			xedges = np.arange(j0, j1, (j1 - j0) / 21.0)
			yedges = np.arange(i0, i1, (i1 - i0) / 21.0)

			H, xedges, yedges = np.histogram2d(pj, pi, bins=(xedges, yedges))

			extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]]
			
			HT1 = []

			HTot = 0.0
			for ii in range(len(xedges)-1):
				for jj in range(len(yedges)-1):
					HT1.append(H[ii][jj])
					HTot += H[ii][jj]


			HT2 = np.sort(HT1)

			nH = len(HT2)
			Hs = 0.0

			HL = []
			for ii in range(nH - 1, -1, -1):
				HOld = Hs
				Hs += HT2[ii]
				if(HOld /HTot < 0.68 and Hs / HTot >= 0.68):
					HL.append(HT2[ii])
				if(HOld /HTot < 0.95 and Hs / HTot >= 0.95):
					HL.append(HT2[ii])
				if(HOld /HTot < 0.997 and Hs / HTot >= 0.997):
					HL.append(HT2[ii])


			HHL = []
			HHL.append(HL[0])
			HHL.append(HL[0]*1.01)
			cs = pl.contour(H.T, colors = 'b', linewidths = 0.6, extent=extent, levels = HHL)
	
			HHL = []
			HHL.append(HL[1])
			HHL.append(HL[1]*1.01)
			cs = pl.contour(H.T, colors = 'b', linewidths = 0.6, extent=extent, levels = HHL)

			#HHL = []
			#HHL.append(HL[2])
			#HHL.append(HL[2]*1.01)
			#cs = pl.contour(H.T, colors = 'b', linewidths = 0.6, extent=extent, levels = HHL)

			
			
			#pl.scatter(pj, pi, s = 0.05, c = col, edgecolors='none')
			print("x", i, xl, j0, j1)
			print("y", j, yl, i0, i1)
		else:
			pl.hist(pi, 21, normed=1, edgecolor = 'none', facecolor='green', alpha=0.75)
			
			pm = np.median(pi)
			pn = np.percentile(pi, 100.0 - 84.1345);
			pp = np.percentile(pi, 84.1345);
			pl.axvline(pn, color='k', linewidth = 0.3)
			pl.axvline(pm, color='k', linestyle='dashed')
			pl.axvline(pp, color='k', linewidth = 0.3)
		
			print("x", i, xl, j0, j1, pm, pm - pn, pp - pm)
			print("y", j, yl, i0, i1)

		if(i == n - 1):
			pl.xlabel(xl)
		if(j == 0):
			pl.ylabel(yl)

		
		pl.xlim([j0,j1])
		if(i != j):
			pl.ylim([i0,i1])

		pl.xticks([])
		pl.yticks([])


name = '../corner_ttv.png'
pl.savefig(name, format='png', dpi=150,bbox_inches='tight')
