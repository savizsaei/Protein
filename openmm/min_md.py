from pdbfixer import PDBFixer
from openmm.app import (PDBFile, Modeller, ForceField, PME, HBonds,
                        DCDReporter, StateDataReporter, Simulation)
from openmm import LangevinMiddleIntegrator, Platform, unit

fixer = PDBFixer(pdbid='1CRN')
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingHydrogens(pH=7.0)

with open("1CRN_fixed.pdb", "w") as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f)

pdb = PDBFile("1CRN_fixed.pdb")
ff = ForceField('amber14-all.xml', 'amber14/tip3p.xml')
mod = Modeller(pdb.topology, pdb.positions)
mod.addSolvent(ff, model='tip3p', padding=1.0*unit.nanometer, ionicStrength=0.05*unit.molar)

# >>> Save solvated topology that matches the trajectory <<<
with open("solvated.pdb", "w") as f:
    PDBFile.writeFile(mod.topology, mod.positions, f)

system = ff.createSystem(mod.topology, nonbondedMethod=PME,
                         nonbondedCutoff=1.0*unit.nanometer, constraints=HBonds)
integrator = LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.004*unit.picoseconds)
platform = Platform.getPlatformByName('CPU')

sim = Simulation(mod.topology, system, integrator, platform)
sim.context.setPositions(mod.positions)

sim.minimizeEnergy()
sim.reporters.append(StateDataReporter("openmm_log.csv", 100, step=True,
    potentialEnergy=True, temperature=True, density=True, speed=True, progress=True,
    totalSteps=5000, separator=','))
sim.reporters.append(DCDReporter("traj.dcd", 100))
sim.step(5000)

print("Done. Wrote solvated.pdb, traj.dcd, openmm_log.csv, 1CRN_fixed.pdb")