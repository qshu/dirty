import rebound
sim = rebound.Simulation()
#set unit
sim.units = ('s', 'au', 'Msun')
# pc = 3.0857 * 10**13 km
# pc = 2.0626 * 10**5 AU
# AU = 1.4960 * 10**8 km
pc2au = 2.0626 * 10**5
au2km = 1.4960 * 10**8

n = 10
leastA = 100.
Comets = []

print 'try to add ',n,'comets ...'
#read Comets.txt
fileName = 'Stars.txt'
File = open(fileName,'r')
lines = File.readlines()[1:]
File.close()

Stars = []
#Mass_[Msun]    x_[pc]  y_[pc]  z_[pc]  vx_[km/s]  vy_[km/s]  vz_[km/s]
for line in lines:
    strs = line.strip().split()
    m = float(strs[0])
    x = float(strs[1]) * pc2au
    y = float(strs[2]) * pc2au
    z = float(strs[3]) * pc2au
    vx = float(strs[4]) / au2km
    vy = float(strs[5]) / au2km
    vz = float(strs[6]) / au2km
    Stars.append([m, x, y, z, vx, vy, vz])

for p in Stars:
    [m, x, y, z, vx, vy, vz] = p 
    sim.add(m=m, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz)
    for i in range(1,n+1):
        sim.add(primary=sim.particles[0], m=10**-10., a= leastA * (1.+ 10*float(i)/n))
        comet = sim.particles[i]
        Comets.append([comet.m, comet.x, comet.y, comet.z, comet.vx, comet.vy, comet.vz])
    #sim.status()
    sim.move_to_com()
    fig = rebound.OrbitPlot(sim, unitlabel="[AU]", color=True)
    fig.savefig('../analysis/Comets.pdf')
    break

# write out to data
fileName = 'data'
File = open(fileName,'w')
for p in Stars:
    File.write(str(p[0]) + ' ' + str(p[1] / pc2au) + ' ' + str(p[2] / pc2au) + ' ' + str(p[3] / pc2au) + ' ' + str(p[4] * au2km) + ' ' + str(p[5] * au2km) + ' ' + str(p[6] * au2km) + '\n') 
for p in Comets:
    File.write(str(p[0]) + ' ' + str(p[1] / pc2au) + ' ' + str(p[2] / pc2au) + ' ' + str(p[3] / pc2au) + ' ' + str(p[4] * au2km) + ' ' + str(p[5] * au2km) + ' ' + str(p[6] * au2km) + '\n')
    #File.write(str(p[0]) + ' ' + str(p[1]) + ' ' + str(p[2]) + ' ' + str(p[3]) + ' ' + str(p[4]) + ' ' + str(p[5]) + ' ' + str(p[6]) + '\n')
File.close()
