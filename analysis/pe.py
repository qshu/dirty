
file = open('../dirty/log.run','r')
runLog = file.readlines()

lines = []
for i in range(len(runLog)):
    line = runLog[i]
    if line[:4] == '  PE':
        peLine = line
        lines.append(runLog[i+1])
def TansferToFloat(i):
    try:
        return float(i)
    except:
        return float(i.replace('D','E'))

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.size'] = 9.0

plt.figure(figsize=(16,9))
labels = peLine.strip().split()[2:]
labels = labels[3:18]

T = []
#print peLine
for line in lines:
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
plt.savefig('PE.pdf')
