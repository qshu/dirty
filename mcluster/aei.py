import sys
rvir = float(sys.argv[1])
tscale = float(sys.argv[2])
NAME = sys.argv[3]

import rebound
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
    if strs[0] == 'N':
        N = int(strs[2])
    elif strs[0] == 'nbin':
        nbin = int(strs[2])
    elif strs[0] == 'units':
        units = int(strs[2])



def readLine(line):
    strs = line.strip().split()
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


orbitInfo = []
File = open(NAME+'.dat.10','r')
for i in range(nbin):
    line = File.readline()
    star = readLine(line)
    line = File.readline()
    comet = readLine(line)

    sim = rebound.Simulation()
    #set unit
    sim.units = ('s', 'au', 'Msun')
    [m, x, y, z, vx, vy, vz] = star
    sim.add(m=m, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz)
    [m, x, y, z, vx, vy, vz] = comet
    sim.add(m=m, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz)
    #sim.status()
    #print sim.particles[1].calculate_orbit(primary=sim.particles[0])
    orbits = sim.calculate_orbits()
    for o in orbits:
        orbitInfo.append([o.a, o.e, o.inc, o.Omega, o.omega, o.f])

    #break
    #sim.move_to_com()
    #fig = rebound.OrbitPlot(sim, unitlabel="[AU]", color=True)
    #fig.savefig('Comets.pdf')
File.close()

orbitInfo = np.array(orbitInfo)

a = orbitInfo[:,0]
e = orbitInfo[:,1]


fig = plt.figure(figsize=(16,9))
ax = fig.add_subplot(111)
ax.semilogx(a,e,'*')
plt.xlabel('semi-major axis (AU)')
plt.ylabel('eccentricity')
plt.title('a vs e')
plt.savefig('../analysis/download/init_loga_e.pdf')

fig = plt.figure(figsize=(16,9))
ax = fig.add_subplot(111)
ax.plot(a,e,'r*')
plt.xlabel('semi-major axis (AU)')
plt.ylabel('eccentricity')
plt.title('a vs e')
plt.savefig('../analysis/download/init_a_e.pdf')
exit()



# write out to data
fileName = 'data'
File = open(fileName,'w')
for p in Stars:
    File.write(str(p[0]) + ' ' + str(p[1] / pc2au) + ' ' + str(p[2] / pc2au) + ' ' + str(p[3] / pc2au) + ' ' + str(p[4] * au2km) + ' ' + str(p[5] * au2km) + ' ' + str(p[6] * au2km) + '\n') 
for p in Comets:
    File.write(str(p[0]) + ' ' + str(p[1] / pc2au) + ' ' + str(p[2] / pc2au) + ' ' + str(p[3] / pc2au) + ' ' + str(p[4] * au2km) + ' ' + str(p[5] * au2km) + ' ' + str(p[6] * au2km) + '\n')
    #File.write(str(p[0]) + ' ' + str(p[1]) + ' ' + str(p[2]) + ' ' + str(p[3]) + ' ' + str(p[4]) + ' ' + str(p[5]) + ' ' + str(p[6]) + '\n')
File.close()
