from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, DataStructs
from rdkit.Chem import rdMolDescriptors as rdmd
from rdkit.Chem import SDWriter

smiles = [
    ("aspirin","CC(=O)OC1=CC=CC=C1C(=O)O"),
    ("caffeine","Cn1cnc2n(C)c(=O)n(C)c(=O)c12"),
    ("ibuprofen","CC(C)CC1=CC=C(C=C1)C(C)C(=O)O")
]
mols = []
for name, smi in smiles:
    m = Chem.MolFromSmiles(smi)
    m = Chem.AddHs(m)
    AllChem.EmbedMolecule(m, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(m)
    m.SetProp("_Name", name)
    mols.append(m)

# fingerprints + similarity
fps = [rdmd.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) for m in mols]
for i in range(len(mols)):
    row = []
    for j in range(len(mols)):
        row.append(DataStructs.TanimotoSimilarity(fps[i], fps[j]))
    print(mols[i].GetProp("_Name"), row)

# basic descriptors
for m in mols:
    print(m.GetProp("_Name"),
          "MW=", Descriptors.MolWt(m),
          "LogP=", Descriptors.MolLogP(m),
          "HBD=", rdmd.CalcNumHBD(m),
          "HBA=", rdmd.CalcNumHBA(m))

# write SDF for later steps
w = SDWriter("ligands.sdf")
for m in mols: w.write(m)
w.close()
print("Wrote ligands.sdf")
