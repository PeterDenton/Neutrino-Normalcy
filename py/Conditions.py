import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage
import nufit
import os

# make output directory
if not os.path.exists("fig"):
    os.makedirs("fig")

# get nufit contours
nf_s23sqs, nf_deltas, nf_chisqs = nufit.get_chisq()
nf_t23s = np.arcsin(np.sqrt(nf_s23sqs)) * 180 / np.pi
# make cosd 1D array
nf_cosds = []
for i in xrange(len(nf_deltas)):
	cosd = np.cos(nf_deltas[i] * np.pi / 180)
	if nf_deltas[i] > 180: continue
	nf_cosds.append(cosd)
nf_cosds = np.array(nf_cosds)
# make chisq 2D array for cosd
nf_chisqs_cosd = np.empty((len(nf_s23sqs), len(nf_cosds)))
for i in xrange(len(nf_t23s)):
	for j in xrange(len(nf_cosds)):
		nf_chisqs_cosd[i, j] = min(nf_chisqs[i, j], nf_chisqs[i, -1 - j])

plt.rc("text.latex", preamble = r"\usepackage{amsmath}")

colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

# nufit 4.1, july 2019
t12 =33.82 * np.pi / 180
t13 =8.61 * np.pi / 180
d = 222. * np.pi / 180
t23 = 48.3 * np.pi / 180

# cupdate the absolute values of the elements of the lepton mixing matrix
def recalc():
	global c12, s12, c13, s13, c23, s23, eid
	global Ue1, Ue2, Ue3, Um1, Um2, Um3, Ut1, Ut2, Ut3
	c12=np.cos(t12)
	s12=np.sin(t12)
	c13=np.cos(t13)
	s13=np.sin(t13)
	c23=np.cos(t23)
	s23=np.sin(t23)
	eid=np.exp(1j*d)
	Ue1=c12*c13
	Ue2=s12*c13
	Ue3=s13
	Um1=abs(-s12*c23-c12*s23*s13*eid)
	Um2=abs(c12*c23-s12*s23*s13*eid)
	Um3=s23*c13
	Ut1=abs(s12*s23-c12*c23*s13*eid)
	Ut2=abs(-c12*s23-s12*c23*s13*eid)
	Ut3=c23*c13

recalc()

# size of the array
n = int(1e3)

# the 1d arrays
t23s = np.linspace(30, 60, n)
s23sqs = np.linspace(0.25, 0.75, n)
deltas = np.linspace(0, 360, n)
cosds = np.linspace(-1, 1, n)

# plot everything
# if t23_b=True then the axis is t23, else it is s23sq
# if cosd_b=True then the yaxis is cosd, else it is just delta
def p(t23_b, cosd_b):
	global t23, d
	cond_taus = np.zeros((n, n))
	cond_mus = np.zeros((n, n))
	cond_1s = np.zeros((n, n))
	cond_2s = np.zeros((n, n))
	cond_3s = np.zeros((n, n))

	for i in xrange(n): # t23
		if t23_b:	t23 = t23s[i] * np.pi / 180
		else:		t23 = np.arcsin(np.sqrt(s23sqs[i]))
		for j in xrange(n): # delta
			if cosd_b:	d = np.arccos(cosds[j])
			else:		d = deltas[j] * np.pi / 180
			recalc()

			# check the conditions
			if Um1 < Um2 and Um3 < Um2:	cond_mus[j, i] = 1.
			if Ut1 < Ut2 and Ut2 < Ut3:	cond_taus[j, i] = 1.
			if Ue1 > Um1 and Um1 > Ut1: cond_1s[j, i] = 1.
			if Um2 > Ue2 and Um2 > Ut2: cond_2s[j, i] = 1.
			if Ue3 < Um3 and Um3 < Ut3: cond_3s[j, i] = 1.

			# if all five normalcy conditions are met and theta23 is "large" (close to the maximal), print theta23 and delta
#			if cond_mus[j, i] + cond_taus[j, i] + cond_1s[j, i] + cond_2s[j, i] + cond_3s[j, i] == 5 and t23 > 0.71: print t23 * 180 / np.pi, d * 180 / np.pi

	if t23_b:
		if cosd_b:	Xs, Ys = np.meshgrid(t23s, cosds)
		else:		Xs, Ys = np.meshgrid(t23s, deltas)
	else:
		if cosd_b:	Xs, Ys = np.meshgrid(s23sqs, cosds)
		else:		Xs, Ys = np.meshgrid(s23sqs, deltas)

	levels = [0.5, 1.5]
	alpha = 0.2
	cs = 1
	lw = 0.5
	z = -20

	# put the labels in the right places

	# condition mu
	plt.contourf(Xs, Ys, cond_mus, levels = levels, alpha = alpha, colors = colors[cs], zorder = z)
	plt.contour(Xs, Ys, cond_mus, levels = levels, colors = colors[cs], linewidths = lw, zorder = z)
	if t23_b:	x = 39.2
	else:		x = 0.4
	if cosd_b:	y = 0.4
	else:		y = 290
	plt.text(x, y, r"$\boldsymbol\mu$", ha = "center", va = "center", zorder = 10)

	# condition tau
	plt.contourf(Xs, Ys, cond_taus, levels = levels, alpha = alpha, colors = colors[cs], zorder = z)
	plt.contour(Xs, Ys, cond_taus, levels = levels, colors = colors[cs], linewidths = lw, zorder = z)
	if t23_b:
		if cosd_b:
			x = 48.1
			y = 0.7
		else:
			x = 52.1
			y = 180
	else:
		if cosd_b:
			x = 0.585
			y = 0
		else:
			x = 0.625
			y = 180
	plt.text(x, y, r"$\boldsymbol\tau$", ha = "center", va = "center", zorder = 10)

	# condition 1
	plt.contourf(Xs, Ys, cond_1s, levels = levels, alpha = alpha, colors = colors[cs], zorder = z)
	plt.contour(Xs, Ys, cond_1s, levels = levels, colors = colors[cs], linewidths = lw, zorder = z)
	if t23_b:	x = 53
	else:		x = 0.64
	if cosd_b:	y = 0.6
	else:		y = 304
	plt.text(x, y, r"{\bf 1}", ha = "center", va = "center", zorder = 10)

	# condition 2
	plt.contourf(Xs, Ys, cond_2s, levels = levels, alpha = alpha, colors = colors[cs], zorder = z)
	plt.contour(Xs, Ys, cond_2s, levels = levels, colors = colors[cs], linewidths = lw, zorder = z)
	if t23_b:	x = 40
	else:		x = 0.41
	if cosd_b:	y = 0.87
	else:		y = 25
	plt.text(x, y, r"{\bf 2}", ha = "center", va = "center", zorder = 10)

	# condition 3
	plt.contourf(Xs, Ys, cond_3s, levels = levels, alpha = alpha, colors = colors[cs], zorder = z)
	plt.contour(Xs, Ys, cond_3s, levels = levels, colors = colors[cs], linewidths = lw, zorder = z)
	if t23_b:	x = 45
	else:		x = 0.5
	if cosd_b:	y = -0.5
	else:		y = 150
	plt.text(x, y, r"{\bf 3}", ha = "center", va = "center", zorder = 10)

	s = r'${\rm Shaded}\Rightarrow$'
	s += r'${\rm\ normal}$'
	if t23_b:	x = 30.2
	else:		x = s23sqs[0] + 0.002
	if cosd_b:	y = 0.98
	else:		y = 356
	plt.text(x, y, s, ha = "left", va = "top", fontsize = 12)

	s = r"${\rm Denton}$"
	s += "\n"
	s += r"$2020$"
	if t23_b:	x = 59.9
	else:		x = s23sqs[-1] - 0.002
	if cosd_b:	y = -0.99
	else:		y = 0.1
	plt.text(x, y, s, ha = "right", va = "bottom", fontsize = 10, color = "darkgray")

	# nufit best fit point
	t23bf = 48.3
	dbf = 222
	nf_c = "gray"
	nf_c = (0.3, 0.3, 0.3)
	nf_c = "b"
	nf_cs = [nf_c]
	if t23_b:	x = t23bf
	else:		x = np.sin(t23bf * np.pi / 180) ** 2
	if cosd_b:	y = np.cos(dbf * np.pi / 180)
	else:		y = dbf
	plt.plot([x], [y], c = nf_c, marker = "+")

	# draw nufit contours
	if t23_b:
		nf_xs = nf_t23s
		x1 = 48
		x2 = 48.5
	else:
		nf_xs = nf_s23sqs
		x1 = 0.553
		x2 = 0.56
	if cosd_b:
		nf_ys = nf_cosds
		y1 = 0.34
		y2 = 0.91
	else:
		nf_ys = nf_deltas
		y1 = 291
		y2 = 348
	nf_Xs, nf_Ys = np.meshgrid(nf_xs, nf_ys)
	s1 = r"$1\sigma$"
	s2 = r"$2\sigma$"
	if cosd_b:	nf_Zs = nf_chisqs_cosd
	else:		nf_Zs = nf_chisqs
	levels = [2.3, 6.18]
	plt.contour(nf_Xs, nf_Ys, nf_Zs.transpose(), levels = levels, linestyles = ":", colors = nf_cs, linewidths = 0.8)
	plt.text(x1, y1, s1, color = nf_c, fontsize = 8, ha = "center", va = "bottom")
	plt.text(x2, y2, s2, color = nf_c, fontsize = 8, ha = "center", va = "bottom")

	# ticks
	if cosd_b:
		plt.yticks([-1, 0, 1])
		plt.gca().set_yticks([-1, -0.5, 0, 0.5, 1], minor = True)
	else:
		plt.yticks([0, 90, 180, 270, 360])

	# label axes
	if t23_b:	plt.xlabel(r"$\theta_{23}{\rm\ [^\circ]}$")
	else:		plt.xlabel(r"$s_{23}^2$")
	if cosd_b:	plt.ylabel(r"$\cos\delta$", labelpad = -10)
	else:		plt.ylabel(r"$\delta{\rm\ [^\circ]}$")

	# make sure axes are limited correctly
	v = [0, 0, 0, 0]
	if t23_b:
		v[0] = t23s[0]
		v[1] = t23s[-1]
	else:
		v[0] = s23sqs[0]
		v[1] = s23sqs[-1]
	if cosd_b:
		v[2] = cosds[0]
		v[3] = cosds[-1]
	else:
		v[2] = deltas[0]
		v[3] = deltas[-1]
	plt.axis(v)

	# filename
	s = "fig/Conditions_"
	if t23_b:	s += "t23_"
	else:		s += "s23sq_"
	if cosd_b:	s += "cosd"
	else:		s += "d"
	plt.savefig(s + ".pdf")

	# analytic approximations for s23sq, cosd graph only
	if not t23_b and cosd_b:
		c = "orange"
		ls = "-"
		lw = 0.5
		xs = np.linspace(0.25, 0.75, 100)
		z = -10

		ys = ((2 * s12 ** 2) * xs - s12 ** 2) / (4 * c12 * s12 * s13 * 0.5) # for condition 1
		plt.plot(xs, ys, c = c, ls = ls, lw = lw, zorder = z)

		ys = ((c12 ** 2 - 0*s12 ** 2 * s13 ** 2) * (1 - 2 * xs)) / (4 * c12 * s12 * s13 * 0.5) # for condition 2
		plt.plot(xs, ys, c = c, ls = ls, lw = lw, zorder = z)

		plt.plot([0.5, 0.5], [-1, 1], c = c, ls = ls, lw = lw, zorder = z) # for condition 3 (exact)

		ys = ((s12 ** 2 * s13 ** 2 - c12 ** 2 - c13 ** 2) * xs + c12 ** 2) / (2 * c12 * s12 * s13 * 0.5) # for condition mu
		plt.plot(xs, ys, c = c, ls = ls, lw = lw, zorder = z)

		ys = ((c13 ** 2 - s12 ** 2 * s13 ** 2 + c12 ** 2) * xs + s12 ** 2 * s13 ** 2 - c13 ** 2) / (-2 * c12 * s12 * s13 * 0.5) # for condition tau
		plt.plot(xs, ys, c = c, ls = ls, lw = lw, zorder = z)

		plt.savefig(s + "_anal.pdf")

	plt.clf()

p(True, True) # t23, cosd
p(False, True) # s23sq, cosd
p(True, False) # t23, delta
p(False, False) # s23sq, delta

