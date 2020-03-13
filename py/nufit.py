import numpy as np

# returns s23sqs, deltas (deg), chisqs[s23sq,delta]
def get_chisq():
	fname = "input/v41.release-SKoff-NO.txt"
	dataf = open(fname, "r")

	s23sqs = []
	deltas = []

	# size of the array
	ns23sq = 101
	ndelta = 109 # -180->360
	ndelta = 73 # 0->360

	chisqs = np.empty((ns23sq, ndelta))

	reading = False
	count = 0
	for line in dataf.readlines():
		if line == "\n": continue
		if line[0] == "#":
			reading = False
			line = line.split()

			# we're only interested in the t23/delta projection
			if line[1] == "T23/DCP":
				reading = True
			continue

		if reading:
			count += 1
			line = line.split()
			s23sq = float(line[0])
			delta = float(line[1])
			chisq = float(line[2])
			if delta < 0: continue # we ignore the delta<0 cases
			if s23sq not in s23sqs: s23sqs.append(s23sq)
			if delta not in deltas: deltas.append(delta)
			is23sq = s23sqs.index(s23sq)
			idelta = deltas.index(delta)
			chisqs[is23sq, idelta] = chisq
	print "nufit read in:"
	print "ns23sq =", len(s23sqs)
	print "ndelta =", len(deltas)

	return np.array(s23sqs), np.array(deltas), chisqs

