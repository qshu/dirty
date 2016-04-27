import h5py
import glob
out = h5py.File("snap.40_total.h5part", "w")
i = 0
shortName = "Step#"
h5files = glob.glob('../dirty/*.h5part')
NBtimes = len(h5files)-1
#print NBtimes,'files'
for t in range(NBtimes):
    fs = h5py.File('../dirty/snap.40_'+str(t)+'.h5part','r')
    steps = len(fs)
    #print steps,'steps'
    #print fs.keys()
    for s in range(steps):
        groupName = shortName + str(s)
        #print groupName
        group_path = fs[groupName].parent.name
        group_id = out.require_group(group_path)
        #print group_id
        newName = shortName + str(i)
        #print newName
        fs.copy(groupName, group_id,name=newName)
        #fs.copy(name, group_id,name=n)
        i = i + 1
    fs.close()
#print out.keys()
out.close()
