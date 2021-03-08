import sys 
import numpy as np

from rdkit import Chem

fp_case = 0
fp_len = 128

  
filename = str(sys.argv[1])
print("read SMILES from file:", filename)
smiles = np.loadtxt(filename, dtype=str)


from rdkit.Chem.AtomPairs.Sheridan import GetBPFingerprint
from rdkit.Chem.EState.Fingerprinter import FingerprintMol
from rdkit.Avalon.pyAvalonTools import GetAvalonFP
from rdkit.Chem.AllChem import  GetMorganFingerprintAsBitVect, GetErGFingerprint
from rdkit.DataStructs.cDataStructs import ConvertToNumpyArray

def generate_fingerprint(mol, fp_case, fp_len):
    print ("Using the fingerprint:")
    if (fp_case==0):
        print ("Estate (1995)")
        return FingerprintMol(mol)[0]
    elif (fp_case==1):
        print ("Morgan circular")
        return GetMorganFingerprintAsBitVect(mol, 2, fp_len)
    elif (fp_case==2):
        print ("Atom pair (1985)")
        return GetHashedAtomPairFingerprintAsBitVect(mol, fp_len)
    elif (fp_case==3):
        print ("Topological torsion (1987)")
        return GetHashedTopologicalTorsionFingerprintAsBitVect(mol, fp_len)
    elif (fp_case==4):
        print ("Avalon bit based (2006)")
        return GetAvalonFP(mol, fp_len)
    elif (fp_case==5):
        print ("Avalon+mol. weight")
        return np.append(GetAvalonFP(mol, fp_len), Descriptors.MolWt(mol))
    elif (fp_case==6):
        print ("RDKit fingerprint")
        return RDKFingerprint(mol, fpSize=fp_len)
    elif (fp_case==7):
        print ("ErG fingerprint (2006)")
        return GetErGFingerprint(mol)


import glob
import pickle

for tsmiles in smiles:
    print ("**********************************")
    print (tsmiles)
    tmol = Chem.MolFromSmiles(tsmiles)
    fp = generate_fingerprint(tmol, fp_case, fp_len)
    nlen = len(fp)
    
    X = np.array(list(fp))
    X = np.reshape(X, (1, nlen)) 
    print ("The shape of the fingerprint is", X.shape)
    
    postname = "_" + str(fp_case) + "_" + str(fp_len) + ".sav"
    files_model = glob.glob('*'+postname)
    for f in files_model:
        model = pickle.load(open(f, 'rb'))
        y_pred = model.predict(X)
        print (f, y_pred)




