print 'transfer to nbody unit ...'
fileName = 'data'
f = open(fileName,'r')
lines = f.readlines()
f.close()

List = []
for line in lines:
	strs = line.strip().split()
	List.append([float(strs[0]), float(strs[1]), float(strs[2]), float(strs[3]), float(strs[4]), float(strs[5]), float(strs[6])])

def DistanceOF(p1,p2):
	return ( (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2 )**0.5

def output(INPUT):
	[M, T, U, W, V, X1, X2, X3, V1, V2, V3] = INPUT
	E = T + U + W
	print '++++++++++++++++++++++++++++++++++'
	print 'Particle numbers:', len(List)
	print 'M:',M
	print 'T:',T
	print 'U:',U
	print 'W:',W
	print 'total energy E = ',E
	print 'virial energy V = ',V
	#print 'Qv', (Qvir*abs(V)/T)**0.5
	#print 'beta:',beta
	print 'X1:',X1
	print 'X2:',X2
	print 'X3:',X3
	print 'V1:',V1
	print 'V2:',V2
	print 'V3:',V3
	print '++++++++++++++++++++++++++++++++++'
	print ''

def calculate(List):
	M = 0. # total mass
	T = 0. # kinetic energy
	U = 0. # potential
	W = 0. # tidal energy
	V = 0. # virial energy


	X1 = 0.
	X2 = 0.
	X3 = 0.
	V1 = 0.
	V2 = 0.
	V3 = 0.

	PL = []
	for i in range(len(List)):
		L = List[i]
		[m,x1,x2,x3,v1,v2,v3] = L

		M += m
		X1 += m*x1
		X2 += m*x2
		X3 += m*x3
		V1 += m*v1
		V2 += m*v2
		V3 += m*v3
		T = T + 0.5*m*(v1**2+v2**2+v3**2)
		for p in PL:
			U += -p[0]*m/DistanceOF(p[1:4],[x1,x2,x3])
		PL.append(L)
	X1 = X1/M
	X2 = X2/M
	X3 = X3/M
	V1 = V1/M
	V2 = V2/M
	V3 = V3/M
	V = U + 2*W

	return([M, T, U, W, V, X1, X2, X3, V1, V2, V3])


data = calculate(List)
print 'before scaling:'
output(data)

[M, T, U, W, V, X1, X2, X3, V1, V2, V3] = data


for i in range(len(List)):
	List[i][1] = List[i][1] - X1
	List[i][2] = List[i][2] - X2
	List[i][3] = List[i][3] - X3
	List[i][4] = List[i][4] - V1
	List[i][5] = List[i][5] - V2
	List[i][6] = List[i][6] - V3

for i in range(len(List)):
	List[i][0] = List[i][0]/M

T = 0. # kinetic energy
U = 0. # potential
W = 0. # tidal energy
V = 0. # virial energy


PL = []
for i in range(len(List)):
	L = List[i]
	[m,x1,x2,x3,v1,v2,v3] = L
	T = T + 0.5*m*(v1**2+v2**2+v3**2)
	for p in PL:
		U += -p[0]*m/DistanceOF(p[1:4],[x1,x2,x3])
	PL.append(L)
V = U + 2*W

Qvir = 0.5 # specified virial ratio
E0 = -0.25 # Nbody unit
Qv = (Qvir*abs(V)/T)**0.5
beta = (1 - Qvir)*U/E0
print 'scaling...'
print 'Qvir = ',Qvir
print 'Qv = ',Qv
print 'beta = ',beta
for i in range(len(List)):
	List[i][4] = List[i][4] * Qv
	List[i][5] = List[i][5] * Qv
	List[i][6] = List[i][6] * Qv

for i in range(len(List)):
	List[i][1] = List[i][1] * beta
	List[i][2] = List[i][2] * beta
	List[i][3] = List[i][3] * beta
	List[i][4] = List[i][4] / (beta)**0.5
	List[i][5] = List[i][5] / (beta)**0.5
	List[i][6] = List[i][6] / (beta)**0.5


data = calculate(List)
print 'after scaling:'
output(data)



fileName = '../dirty/dat.10'
f = open(fileName,'w')
#f.write('# NBODY unit\n')
for j in range(len(List)):
	line = List[j]
	for i in range(len(line)):
		s = line[i]
		f.write(str(s))
		if i < len(line) - 1:
			f.write(' ')
	if j < len(List) - 1:
		f.write('\n')
f.close()
