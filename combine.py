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
    #print sort_strings_with_emb_numbers(fs.keys())
    for groupName in sort_strings_with_emb_numbers(fs.keys()):
        #print fs[groupName].keys()
        group_path = fs[groupName].parent.name
        group_id = out.require_group(group_path)
        #print group_id
        newName = shortName + str(i)
        #print newName
        fs.copy(groupName, group_id,name=newName)
        #fs.copy(name, group_id,name=n)
        i = i + 1
    fs.close()
    #break
#print out.keys()
out.close()
