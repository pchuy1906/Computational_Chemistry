import numpy as np
import json

with open('file.json') as f:
  data = json.load(f)

syms =  data['atoms']['elements']['type']
xyzs =  data['atoms']['coords']['3d']


natom = int(len(xyzs)/3)
print ("the number of atoms is ", natom)

xyzs = np.reshape( np.array(xyzs), (natom, 3) )

f = open("output.xyz", "w")
f.write("%-d\n" % (natom) )
f.write("%-s\n" % ("Comment") )
for iatom in range(natom):
    f.write("%s %15.9f %15.9f %15.9f\n" % (syms[iatom], xyzs[iatom,0], xyzs[iatom,1], xyzs[iatom,2] ) )
f.close()
