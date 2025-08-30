# Mini CADD Practice 
**Docking (Smina / Vina‑compatible), short MD (OpenMM), file prep (Open Babel / RDKit), and visualization (PyMOL)**

This README captures a reproducible, end‑to‑end practice workflow we executed:
- **Protein prep** (PDBFixer) → **PDBQT** (Open Babel)
- **Ligand prep** from SMILES (Open Babel / RDKit)
- **Docking** with **Smina** (a Vina fork) and score collection
- **Short solvated MD** with **OpenMM** and **RMSD** analysis (MDTraj)
- **Visualization** with **PyMOL** (cartoon + ligand sticks, pocket, H‑bond guess)

> Tested on macOS with a conda env named `cadd`. GPU not required.

---

## 0) Environment Setup

```bash
# create env (ARM‑friendly); Vina is replaced with Smina here
conda create -n cadd -c conda-forge python=3.10 rdkit openmm openmmforcefields \
  pdbfixer openbabel mdtraj mdanalysis py3Dmol -y
conda activate cadd

# docking engine: Smina (Vina‑compatible)
conda install -c conda-forge smina -y

# (optional) PyMOL open‑source GUI
conda install -c conda-forge pymol-open-source -y
# If GUI fails to start, install XQuartz:  brew install --cask xquartz
```

Directory used below (adjust if you prefer a different path):
```bash
mkdir -p ~/cadd-practice/viva ~/cadd-practice/openmm
```

---

## 1) Receptor Prep (Crambin, PDB: 1CRN → PDB/PDBQT)

```bash
cd ~/cadd-practice/viva
python - << 'PY'
from pdbfixer import PDBFixer
from openmm.app import PDBFile
fixer = PDBFixer(pdbid='1CRN')
fixer.findMissingResidues(); fixer.findMissingAtoms(); fixer.addMissingHydrogens(pH=7.0)
with open("receptor.pdb","w") as f: PDBFile.writeFile(fixer.topology, fixer.positions, f)
PY

# Convert to PDBQT (adds Gasteiger partial charges)
obabel -ipdb receptor.pdb -opdbqt -O receptor.pdbqt --partialcharge gasteiger
```

---

## 2) Ligand Prep (SMILES → 3D SDF → PDBQT)

Create a small SMILES list and generate 3D structures:
```bash
cat > ligands.smi << 'EOF'
CC(=O)OC1=CC=CC=C1C(=O)O aspirin
Cn1cnc2n(C)c(=O)n(C)c(=O)c12 caffeine
CC(C)CC1=CC=C(C=C1)C(C)C(=O)O ibuprofen
EOF

# SMILES → SDF (3D, add H) → PDBQT (Gasteiger)
obabel -ismi ligands.smi -osdf -O ligands.sdf --gen3d -h
obabel -isdf ligands.sdf -opdbqt -O ligands.pdbqt --partialcharge gasteiger

# Also write per‑ligand PDBQTs
obabel -isdf ligands.sdf -opdbqt -O aspirin.pdbqt   -f 1 -l 1 --partialcharge gasteiger
obabel -isdf ligands.sdf -opdbqt -O caffeine.pdbqt  -f 2 -l 2 --partialcharge gasteiger
obabel -isdf ligands.sdf -opdbqt -O ibuprofen.pdbqt -f 3 -l 3 --partialcharge gasteiger
```

(Alternative with **RDKit** for 3D + similarity; optional)
```python
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, DataStructs, rdMolDescriptors as rdmd
smis = {
    'aspirin':  'CC(=O)OC1=CC=CC=C1C(=O)O',
    'caffeine': 'Cn1cnc2n(C)c(=O)n(C)c(=O)c12',
    'ibuprofen':'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O'
}
mols = []
for name, smi in smis.items():
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    AllChem.EmbedMolecule(m, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(m)
    m.SetProp('_Name', name); mols.append(m)
fps = [rdmd.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) for m in mols]
for i,a in enumerate(mols):
    print(a.GetProp('_Name'), [f'{DataStructs.TanimotoSimilarity(fps[i],fps[j]):.2f}' for j in range(len(mols))])
```

---

## 3) Docking with **Smina** (Vina‑compatible)

**3.1** Choose a search‑box center (protein centroid is OK for practice):
```bash
python - << 'PY'
import MDAnalysis as mda
u = mda.Universe('receptor.pdb')
c = u.atoms.positions.mean(axis=0)  # Angstroms
print(f"{c[0]:.3f} {c[1]:.3f} {c[2]:.3f}")
PY
# Suppose it prints: 9.269 9.820 6.864
```

**3.2** Dock one ligand (aspirin) and inspect scores:
```bash
smina --receptor receptor.pdbqt --ligand aspirin.pdbqt \
  --center_x 9.269 --center_y 9.820 --center_z 6.864 \
  --size_x 20 --size_y 20 --size_z 20 \
  --exhaustiveness 8 --out docked_aspirin.pdbqt --log aspirin.log

grep -i 'mode |' -n aspirin.log || true
```

**3.3** Rank three ligands quickly:
```bash
for L in aspirin caffeine ibuprofen; do
  smina --receptor receptor.pdbqt --ligand ${L}.pdbqt \
    --center_x 9.269 --center_y 9.820 --center_z 6.864 \
    --size_x 20 --size_y 20 --size_z 20 \
    --exhaustiveness 8 --num_modes 1 \
    --out docked_${L}.pdbqt --log ${L}.log
  printf "%-10s" "$L"; awk '/^[[:space:]]*1[[:space:]]*-/{print $2,"kcal/mol"}' ${L}.log
done | awk 'BEGIN{print "ligand,affinity_kcal"}1' > docking_scores.csv
cat docking_scores.csv
```
_Example scores observed during practice (seed‑dependent): ibuprofen ~−4.5, aspirin ~−4.0, caffeine ~−3.7 kcal/mol (weak binders as expected for a generic box on crambin)._

**3.4** Convert a pose for viewers (optional):
```bash
obabel -ipdbqt docked_aspirin.pdbqt -osdf -O docked_aspirin.sdf
```

**3.5 (Optional)** Refine around the best pose (autobox to pose):
```bash
smina --receptor receptor.pdbqt --ligand ibuprofen.pdbqt \
  --autobox_ligand docked_ibuprofen.pdbqt --autobox_add 4 \
  --exhaustiveness 16 --num_modes 9 \
  --out docked_ibuprofen_refined.pdbqt --log ibuprofen_refined.log
grep -i affinity ibuprofen_refined.log
```

---

## 4) Short MD with **OpenMM** + RMSD

Create and run a tiny MD (solvated, ~5k steps) and compute RMSD:
```bash
cd ~/cadd-practice/openmm

cat > min_md.py << 'PY'
from pdbfixer import PDBFixer
from openmm.app import (PDBFile, Modeller, ForceField, PME, HBonds,
                        DCDReporter, StateDataReporter, Simulation)
from openmm import LangevinMiddleIntegrator, Platform, unit
fixer = PDBFixer(pdbid='1CRN')
fixer.findMissingResidues(); fixer.findMissingAtoms(); fixer.addMissingHydrogens(pH=7.0)
with open("1CRN_fixed.pdb","w") as f: PDBFile.writeFile(fixer.topology, fixer.positions, f)
pdb = PDBFile("1CRN_fixed.pdb")
ff = ForceField('amber14-all.xml','amber14/tip3p.xml')
mod = Modeller(pdb.topology, pdb.positions)
mod.addSolvent(ff, model='tip3p', padding=1.0*unit.nanometer, ionicStrength=0.05*unit.molar)
with open("solvated.pdb","w") as f: PDBFile.writeFile(mod.topology, mod.positions, f)
system = ff.createSystem(mod.topology, nonbondedMethod=PME,
                         nonbondedCutoff=1.0*unit.nanometer, constraints=HBonds)
integrator = LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.004*unit.picoseconds)
sim = Simulation(mod.topology, system, integrator, Platform.getPlatformByName('CPU'))
sim.context.setPositions(mod.positions)
sim.minimizeEnergy()
sim.reporters.append(StateDataReporter("openmm_log.csv",100,step=True,temperature=True,
    potentialEnergy=True, progress=True, totalSteps=5000, separator=','))
sim.reporters.append(DCDReporter("traj.dcd",100))
sim.step(5000)
print("Done: solvated.pdb, traj.dcd, openmm_log.csv")
PY

python min_md.py

# RMSD on protein atoms only (topology must match trajectory → use solvated.pdb)
python - << 'PY'
import mdtraj as md, numpy as np
t = md.load('traj.dcd', top='solvated.pdb')
tp = t.atom_slice(t.topology.select('protein'))
r = md.rmsd(tp, tp, 0)  # nm
np.savetxt('rmsd.csv', r, delimiter=',')
print(f"Protein RMSD first/median/last (nm): {r[0]:.3f}, {np.median(r):.3f}, {r[-1]:.3f}")
PY
```

_Example we observed: median RMSD ≈ 0.10–0.12 nm over ~20 ps (healthy micro‑drift)._

(Plot, optional)
```python
import numpy as np, matplotlib.pyplot as plt
r = np.loadtxt('rmsd.csv', delimiter=',')
t = np.arange(len(r))*0.4  # ps if DCD stride=100 and dt=0.004 ps
plt.plot(t, r); plt.xlabel('Time (ps)'); plt.ylabel('Protein RMSD (nm)')
plt.tight_layout(); plt.savefig('rmsd.png', dpi=150)
```

---

## 5) Visualize with **PyMOL**

**GUI launch with both structures:**
```bash
cd ~/cadd-practice/viva
pymol receptor.pdb docked_aspirin.sdf   # or docked_aspirin.pdbqt
```

**In PyMOL console:**
```
hide everything
show cartoon, receptor; color slate, receptor
show sticks, docked_aspirin; color yellow, docked_aspirin
zoom docked_aspirin, 8

# (optional) pocket & quick H-bonds (visual diagnostic)
select pocket, byres (receptor within 4 of docked_aspirin)
show sticks, pocket; color salmon, pocket
distance hbonds, (receptor and name N+O), (docked_aspirin and name N+O), 3.2
hide labels, hbonds

# save figure
bg_color white; set ray_opaque_background, off
ray 1600,1200; png pose.png, dpi=200
```

**Headless one‑click render:**
```bash
cat > pose.pml <<'PML'
load receptor.pdb, receptor
load docked_aspirin.sdf, lig
hide everything
show cartoon, receptor; color slate, receptor
show sticks, lig; color yellow, lig
bg_color white; set ray_opaque_background, off
zoom lig, 8
ray 1600,1200
png pose.png, dpi=200
PML
pymol -cq pose.pml
```

---

## What You Can Honestly Claim Now

- **OpenMM** (prep, solvate, minimize, short MD; RMSD analysis)  
- **Open Babel** (format conversions; PDB↔PDBQT; SMILES→3D SDF; Gasteiger charges)  
- **AutoDock Vina–based docking** (via **Smina**) with logs and poses  
- **RDKit** (3D conformers, fingerprints, basic descriptors)  
- **PyMOL** (loading, styling, pocket selection, rendering)

---

## Acknowledgments

- **Smina** (Vina‑compatible docking), **Open Babel**, **OpenMM + PDBFixer**, **RDKit**, **MDTraj**, **PyMOL**.  
Please cite upstream tools in any publication or report.
