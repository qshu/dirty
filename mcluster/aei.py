import sys
rvir = float(sys.argv[1])
tscale = float(sys.argv[2])
NAME = sys.argv[3]
N = int(sys.argv[4])
nbin = int(sys.argv[5])
amin = float(sys.argv[6])
amax = float(sys.argv[7])

import rebound
import random
# pc = 3.0857 * 10**13 km
# pc = 2.0626 * 10**5 AU
# AU = 1.4960 * 10**8 km
pc2au = 2.0626 * 10**5
au2km = 1.4960 * 10**8

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


File = open(NAME+'.info','r')
lines = File.readlines()
File.close()
for line in lines:
    if line == '\n': continue
    strs = line.strip().split()
    if strs[0] == 'units':
        units = int(strs[2])
        #print '............ unit = ', units



def readLine(line):
    strs = line.split()
    m = float(strs[0])
    x = float(strs[1]) * pc2au
    y = float(strs[2]) * pc2au
    z = float(strs[3]) * pc2au
    vx = float(strs[4]) / au2km
    vy = float(strs[5]) / au2km
    vz = float(strs[6]) / au2km
    if units == 0:
        #print 'nbody units --> astrophysical units'
        m *= N
        x *= rvir
        y *= rvir
        z *= rvir
        vx *= rvir/tscale
        vy *= rvir/tscale
        vz *= rvir/tscale
    elif units != 1:
        print 'unknown units'
    return [m, x, y, z, vx, vy, vz]
def writeLine(l):
    re = ''
    for i in l:
        re = re + str(i) + ' '
    return re[:-1] + '\n'


orbitInfo = []
binarylines = []
File = open(NAME+'.dat.10','r')
for i in range(nbin):
    ta = amin + (amax - amin)/float(nbin)* i
    line = File.readline().strip()
    #print line
    star = readLine(line)

    sim = rebound.Simulation()
    sim.units = ('s', 'au', 'Msun')
    [m, x, y, z, vx, vy, vz] = star
    sim.add(m=m, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz)
    #sim.add(m=0, a=1, e=0, inc=0, Omega=0, omega=0, f=0)
    sim.add(m=0., a=ta, e=random.random(), inc=random.random()*np.pi, Omega=random.random()*np.pi, omega=random.random()*np.pi, f=random.random()*np.pi)
    #[m, x, y, z, vx, vy, vz] = comet
    #sim.add(m=m, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz)
    #sim.status()
    #print sim.particles[1].calculate_orbit(primary=sim.particles[0])
    orbit = sim.particles[1].calculate_orbit()
    [xc, yc, zc] = sim.particles[1].xyz
    [vxc, vyc, vzc] = sim.particles[1].vxyz

    sline = [m/N, x/rvir/pc2au, y/rvir/pc2au, z/rvir/pc2au, vx/(rvir/tscale)*au2km, vy/(rvir/tscale)*au2km, vz/(rvir/tscale)*au2km]
    cline = [m/N*10**-10, xc/rvir/pc2au, yc/rvir/pc2au, zc/rvir/pc2au, vxc/(rvir/tscale)*au2km, vyc/(rvir/tscale)*au2km, vzc/(rvir/tscale)*au2km]
    lineStar = writeLine(sline)
    lineComet = writeLine(cline)

    binarylines.append(lineStar)
    binarylines.append(lineComet)


    #print star
    #print orbit
    #print xyz
    #print vxyz
    #print binarylines
    #exit(0)
    #sim.move_to_com()
    #fig = rebound.OrbitPlot(sim, unitlabel="[AU]", color=True)
    #fig.savefig('Comets.pdf')
singleStarlines = File.readlines()

Wlines = binarylines + singleStarlines
File.close()
File = open('dat.10','w')
File.writelines(Wlines)
File.close()

exit(0)
orbitInfo = np.array(orbitInfo)

a = orbitInfo[:,0]
e = orbitInfo[:,1]


fig = plt.figure(figsize=(16,9))
ax = fig.add_subplot(111)
ax.semilogx(a,e,'*')
plt.xlabel('semi-major axis (AU)')
plt.ylabel('eccentricity')
plt.title('a vs e')
plt.savefig('../log/init_loga_e.pdf')

fig = plt.figure(figsize=(16,9))
ax = fig.add_subplot(111)
ax.plot(a,e,'r*')
plt.xlabel('semi-major axis (AU)')
plt.ylabel('eccentricity')
plt.title('a vs e')
plt.savefig('../log/init_a_e.pdf')



# write out to data
fileName = 'data'
File = open(fileName,'w')
for p in Stars:
    File.write(str(p[0]) + ' ' + str(p[1] / pc2au) + ' ' + str(p[2] / pc2au) + ' ' + str(p[3] / pc2au) + ' ' + str(p[4] * au2km) + ' ' + str(p[5] * au2km) + ' ' + str(p[6] * au2km) + '\n') 
for p in Comets:
    File.write(str(p[0]) + ' ' + str(p[1] / pc2au) + ' ' + str(p[2] / pc2au) + ' ' + str(p[3] / pc2au) + ' ' + str(p[4] * au2km) + ' ' + str(p[5] * au2km) + ' ' + str(p[6] * au2km) + '\n')
    #File.write(str(p[0]) + ' ' + str(p[1]) + ' ' + str(p[2]) + ' ' + str(p[3]) + ' ' + str(p[4]) + ' ' + str(p[5]) + ' ' + str(p[6]) + '\n')
File.close()
