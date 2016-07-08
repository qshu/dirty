import h5py
import glob
import re

re_digits = re.compile(r'(\d+)')
def emb_numbers(s):
    pieces=re_digits.split(s)
    pieces[1::2]=map(int,pieces[1::2])
    return pieces

def sort_strings_with_emb_numbers(alist):
    aux = [(emb_numbers(s),s) for s in alist]
    aux.sort()
    return [s for __,s in aux]

out = h5py.File("snap.40_total.h5part", "w")
i = 0
shortName = "Step#"
h5files = glob.glob('./dirty/*.h5part')
NBtimes = len(h5files)-1
print NBtimes,'files'
for t in range(NBtimes):
    fs = h5py.File('./dirty/snap.40_'+str(t)+'.h5part','r')
    print "add file", t,':', sort_strings_with_emb_numbers(fs.keys())
    for groupName in sort_strings_with_emb_numbers(fs.keys()):
        #print fs[groupName].keys()
        group_path = fs[groupName].parent.name
        group_id = out.require_group(group_path)
        #print group_id
        newName = shortName + str(i)
        #print newName
        fs.copy(groupName, group_id,name=newName)


        #print out.keys()
        #fs.copy(name, group_id,name=n)
        i = i + 1
    fs.close()
    #break
print 'out:'
for groupName in sort_strings_with_emb_numbers(out.keys()):
    out[groupName].create_dataset("ID", (len(out[groupName]["X1"]),), 'int32')
    #out[groupName]["ID"] = int(groupName[-1])
    print groupName, out[groupName].keys() , out[groupName]["ID"]
out.close()
#"""
#    Arg-sort the h5part particles with the ID array.
#"""
#
#import h5py
#import numpy as np
#
#h5f = h5py.File('./dirty/snap.40_1.h5part', 'r')
#h5f_out = h5py.File('data_out.h5part', 'w')
#
#nsteps = len(h5f)
#for i in range(nsteps):
#    print 'Processing Step#%d' % i
#    h5s = h5f['/Step#%d' % i]
#    ID = np.array(h5s['NAM'][...], dtype=np.int) - 1
#    X = np.array(h5s['X1'][...])[ID]
#    Y = np.array(h5s['X2'][...])[ID]
#    Z = np.array(h5s['X3'][...])[ID]
#    M = np.array(h5s['M'][...])[ID]
#
#    h5s_out = h5f_out.create_group('/Step#%d' % i)
#    h5s_out.create_dataset('ID', data=ID)
#    h5s_out.create_dataset('X', data=X)
#    h5s_out.create_dataset('Y', data=Y)
#    h5s_out.create_dataset('Z', data=Z)
#    h5s_out.create_dataset('Mass', data=M)
#
#h5f.close()
#h5f_out.close()
