
file = open('../dirty/log.run','r')
runLog = file.readlines()
file.close()

PElines = []
NBlines = []
LGlinesALL = []
LGlinesMASSLESS = []
for i in range(len(runLog)):
    line = runLog[i]
    if line[:4] == '  PE':
        peLine = line
        PElines.append(runLog[i+1])
    elif line[1:9] == 'TIME[NB]':
        NBlines.append(line)
    elif line[11:31] == 'TIME   M/MT(  all  )':
        LGlinesALL.append(runLog[i+1])
    elif line[11:31] == 'TIME  M/MT(massless)':
        LGlinesMASSLESS.append(runLog[i+1])
#print LGlinesALL
#print LGlinesMASSLESS

#exit(0)
print 'reading log finished.'

def TansferToFloat(i):
    try:
        return float(i)
    except:
        return float(i.replace('D','E'))
import matplotlib as mpl
mpl.use('Agg')
#mpl.rcParams['font.size'] = 9.0


import matplotlib.pyplot as plt

import numpy as np
#############################Lagr part begin###############################
'''
TIME   M/MT:  1.00D-03  3.00D-03  5.00D-03  1.00D-02  3.00D-02  5.00D-02  1.00D-01  2.00D-01  3.00D-01  4.00D-01  5.00D-01  6.00D-01  7.00D-01  8.00D-01  9.00D-01  9.50D-01  9.90D-01  1.00D+00  <RC 
'''

def LGoutput(LGlines, string):
    LGinfo = []
    for line in LGlines:
        temp = line.strip().split()
        info = []
        for i in range(len(temp)):
            if i == 1: continue
            info.append(TansferToFloat(temp[i]))
        LGinfo.append(info)
    LGinfo = np.array(LGinfo)
    fig = plt.figure(figsize=(16,9))
    ax = fig.add_subplot(111)
    T = LGinfo[:,0]
    n = np.shape(LGinfo)[1]
    for i in range(1,n):
        ax.semilogy(T, LGinfo[:,i])
    plt.xlabel('Time [NB]')
    plt.ylabel('Lagrangian radii [NB]')
    plt.title(str(n-1)+' Lagrangian radii for '+string+' particles, M/MT from 1e-3 to 1e0')
    plt.savefig('./download/lagr'+string+'.pdf')


LGoutput(LGlinesALL,'ALL')
LGoutput(LGlinesMASSLESS,'MASSLESS')


'''

LGinfoALL = []
LGinfoMASSLESS = []

for line in LGlinesALL:
    temp = line.strip().split()
    info = []
    for i in range(len(temp)):
        if i == 1: continue
        info.append(TansferToFloat(temp[i]))
    LGinfoALL.append(info)


for line in LGlinesMASSLESS:
    temp = line.strip().split()
    info = []
    for i in range(len(temp)):
        if i == 1: continue
        info.append(TansferToFloat(temp[i]))
    LGinfoMASSLESS.append(info)

LGinfoALL = np.array(LGinfoALL)
LGinfoMASSLESS = np.array(LGinfoMASSLESS)



fig = plt.figure(figsize=(16,9))
ax = fig.add_subplot(111)
T = LGinfoALL[:,0]
n = np.shape(LGinfoMASSLESS)[1]
for i in range(1,n):
    ax.semilogy(T, LGinfoALL[:,i])
plt.xlabel('Time [NB]')
plt.ylabel('Lagrangian radii [NB]')
plt.title(str(n-1)+' Lagrangian radii for ALL particles, M/MT from 1e-3 to 1e0')
plt.savefig('./download/lagrALL.pdf')

fig = plt.figure(figsize=(16,9))
ax = fig.add_subplot(111)
T = LGinfoMASSLESS[:,0]
n = np.shape(LGinfoMASSLESS)[1]
for i in range(1,n):
    ax.semilogy(T, LGinfoMASSLESS[:,i])
plt.xlabel('Time [NB]')
plt.ylabel('Lagrangian radii [NB]')
plt.title(str(n-1)+' Lagrangian radii for MASSLESS particles, M/MT from 1e-3 to 1e0')
plt.savefig('./download/lagrMASSLESS.pdf')

'''
#############################Lagr part end#################################



#exit(0)
#############################NB part begin#################################
'''
['TIME[NB]', '0.000000E+00', 'N', '10000', '<NB>', '19', 'NPAIRS', '0', 'NMERGE', '0', 'MULT', '0', 'NS', '9000', 'NSPL=', '1000', 'NSTEP(I,B,R,U)', '0', '0', '0', '0', 'DE', '0.000000E+00', 'E', '-0.249938']

'''
NBinfo = []
for line in NBlines:
    temp = line.strip().split()
    TIME = float(temp[1])
    N = float(temp[3])
    Nbinary = float(temp[5])
    NPAIRS = float(temp[7])
    NMERGE = float(temp[9])
    MULT = float(temp[11])
    NS = float(temp[13])
    NSPL = float(temp[15])
    DE = float(temp[22])
    E = float(temp[24])
    NBinfo.append([TIME, N, Nbinary, NPAIRS, NMERGE, MULT, NS, NSPL, DE, E])
NBinfo = np.array(NBinfo)
#print NBinfo[:,[0,1]]
T = NBinfo[:,0]

Ylabel = ['N', 'Nbinary', 'NPAIRS', 'NMERGE', 'MULT', 'NS', 'NSPL', 'DE', 'E']

for y in range(len(Ylabel)):
    fig = plt.figure(figsize=(16,9))
    ax = fig.add_subplot(111)
    ax.plot(T, NBinfo[:,y+1])
    plt.xlabel('TIME[NB]')
    plt.ylabel(Ylabel[y])
    plt.savefig('./download/'+Ylabel[y]+'.pdf')
#############################NB part end#################################


#############################PE part begin#################################
plt.figure(figsize=(16,9))
labels = peLine.strip().split()[2:]
labels = labels[3:18]

T = []
#print peLine
for line in PElines:
    temp = line.strip().split()[2:]
    times = [TansferToFloat(i) for i in temp]
    times = times[3:18]
    T.append(times)
steps = len(T)

# This first pie chart fill the plot, it's the lowest level
patches, texts= plt.pie(T[steps-1], labels=labels, labeldistance=1.05)
        #colors=plt.style.library['bmh']['axes.color_cycle'])
for text in texts:
    if text._text == 'Total' or text._text == 'Irr.' or text._text == 'Intgrt' or text._text == 'Reg.':
        text.set_fontsize(30)
'''
Total
Inti.
Intgrt

Reg.
Irr.
Pred.
Init.B.
Mdot
Move
Comm.I.
Comm.R.
Send.I.
Send.R.
KS
Adjust
OUT
Barr.

Barr.I.
Barr.R.
Reg.GPU.S
Reg.GPU.P
Comm.Adj.
Mdot.Fic.
Mdot.Fc.
Mdot.Pot.
Mdot.EC.
Sort.B.
HighV
KS.Init.B
KS.Int.S
KS.Int.P
KS.Comm.
KS.Barr.
KS.Move
KS.Cmb.
KS.Insert
KS.Init.
KS.Term.
Hiar.
KS.UP
xtsub1
xtsub2
'''
ax = plt.gca()
# For each successive plot, change the max radius so that they overlay
for i in range(steps-1,0,-1):
    ax.pie(T[i], radius=float(i+1)/float(steps), autopct = '%3.1f%%',
            pctdistance = ((2*(i+1)-1)/float(2*steps))/((i+1)/float(steps)), shadow = False) 
            #colors=plt.style.library['bmh']['axes.color_cycle'])

ax.set_aspect('equal')

#plt.legend(labels, loc=(1.3, -0.1), shadow=True)
plt.legend(labels, loc=(1.2, 0.3), shadow=True)
#plt.show()
plt.savefig('./download/PE.pdf')

#############################PE part end#################################
