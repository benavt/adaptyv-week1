#!/usr/bin/env python

from pdbfixer import PDBFixer
from openmm import *
from openmm.app import *
from openmm import unit
import os
import re 
import time
import datetime
import numpy as np
import math 

"""
  A program to minimize a PDB file/s with OpenMM (RTT, March 2023. V1.0)

  Args:
       -af2                 Make all defaults to emulate AF2 relax 
       -afw2                Same as above with implicit water 
       -rtt                 Deafults to RT method

       -inD <string>        directory to get the PDB files from
       -outD <string>       directory write final minimized PDB files 
       -list <string>       read from List the PDB file names to minimize
       -steps <int>         the number of steps for minimization,    default 1000
       -tolerance <float>   tolerance in forces for minimization,    default 10.0
       -Ecut <float>        Epot cutoff to report a list of good E ones
       -ffield <string>     Force field to use,                      default AMBER14
       -water <string>      Model for water molecules,               default AMBER14
       -box <float>         Size of cubic box of explicit water (nm),    def: 1.43  
       -restr <string>      restrain atoms to original positions,    
       -cons <string>       constrain atoms (not move at all),       
       -kf  <float>         value of k for force to restrain atoms,  default 100.0
       -NBmethod <string>   nonbonded method to use.     default: CutoffNonPeriodic 
                            choices: NoCutoff, CutoffNonPeriodic, CutoffPeriodic, Ewald, PME
       -NBcut <float>       nonbonded cut-off ")
       -report              Report minimization info/steps
       -gr[oup] <string>    File name for the multimodel file created at the end, if requested

"""

# ----------------------- help message -------------------------------
if ( len(sys.argv) < 2 ):
     print("\n ** OpenMmini.py (RTT, March 2024. V1.0) \n")
     print(" A python script to perform an energy minimization on a PDB file\n")
     print(" Optional Args: options between (), defaults between brackets on the right")
     print("      Special flags                                                           ")
     print("      -af2                   defaults to emulate AF2 relax protocol           ")
     print("      -af2w                  same as af2 plus implicit water in the protocol  ")
     print("      -rtt                   defaults to RTT relax protocol                   ")
     print("                                                                              ")
     print("      -inD <string>          directory to get the PDB files from              ")
     print("      -outD <string>         directory to write final mimized PDB files       ")
     print("      -list <string>         read from List the PDB files to minimize         ")
     print("      -st[eps] <int>         the number of steps for minimization,     [0]    ")
     print("      -tol[erance] <float>   tolerance in forces for minimization,     [10.0] ")
     print("      -Ecut <float>          Epot cutoff to report a list of good E,   [0.0]  ")
     print("      -ff[ield] <string>     Force field to use (amber14, amber99)     [AMBER14] ")
     print("      -wa[ter] <string>      Water model (implicit, amber14, amber99)  [no default]  ")
     print("      -box <float>           Size of cubic box of explicit water (nm), [1.43]        ")
     print("      -res[train] <string>   restrain atoms around original positions, [no default]  ")
     print("      -con[strain] <string>  Constrain atoms, no move at all           [no default]  ")
     print("                             (backbone, heavy, N,CA,C,O  ) ")
     print("      -kf  <float>           value of stiffness to restrain/constrain  [10.0]        ")
     print("      -NBm[ethod] <string>   nonbonded method to use.           [CutoffNonPeriodic]  ")
     print("                             (NoCutoff, CutoffNonPeriodic, CutoffPeriodic) ")
     print("      -NBc[ut] <float>       nonbonded cut-off                         [20 A]")
     print("      -gr[oup] <string>      File name for the multimodel file ")
     print("")
     print(" Defaults for flags: (note than can be changed at command line)           ")
     print("			-af2	    -af2w	 -rtt		no flag       ")
     print("		       --------    -----------  ------------   ------------   ")
     print("	  Method	AF2	    AF2W	RTT		Custom        ")
     print("	  Steps		100	    100		0 (unlimited)	0 (unlimited) ")
     print("	  tolerance	2.39	    2.39	1.0 		1.0           ")
     print("	  restrained	heavy	    heavy	backbone	  --          ")
     print("	  stiffness	10.0	    10.0	5.0     	5.0           ")
     print("	  NBmethod   	NoCutoff    NoCutoff	NoCutoff        NoCutoff      ")
     print("	  ForceField    amber99sb   amber99sb	amber99sb	amber99sb     ")
     print("	  FFwater          --       implicit       --		  --          ")
     print("")
     print(" Use examples:")
     print("      OpenMM_mini.py -af2 -inD Unrelaxed -list ListToMini -outD Minimized -gr AllPDBS.pdb ")
     print("      OpenMM_mini.py -restr heavy -wat implicit  ../SomeDir/*.pdb")
     print("      OpenMM_mini.py -steps 1000 -tolerance 2.39 -list Lista -inD ../SomeDir/")
     print("")
     sys.exit(0)

# ----- General units and conversion constants
ENERGY = unit.kilocalories_per_mole
LENGTH = unit.angstroms
toKjnm2 = 100.0/0.239006
toKjmol = 4.184
toKcalmol = 0.239006

# ----------------------------------- defaults updated at command line
# ----- Special flags, report, plot
AF2 = False 
AFw = False 
RTT = False
Reporting = False   # wether report (or not) on minimization steps
Grouping = False    # wether group (or not) models into a multimodel file at the end
# ----- PDB or mmCIF formats
mmCIF = False 
# ----- steps, tolerance, force field and water ff, integrator
steps = 99      	 # default number of steps, 0 means no limit, changed later
tolerance = 99.0 	 # default in OpenMM, AF2 uses 2.39, changed later
atomsrestrained = ''     # default none
atomsconstrained = ''    # default none
stiffness = 99.0         # Kf for restrained/constrained, changed laterº
Ecut = 1.0e99            # Cutoff for good Epot, changed later
ffield = 'amber99sb.xml'
ffwater = ''
Water = False
# ----- integrators, restrained and force or constrained atoms
integrator = openmm.LangevinIntegrator(0, 0.01, 0.0)
# ----- Others, dir for input/output, cutoff for E 
indirname = './'
outdirname = './'
# -----
files = []
lEp = []
Files_to_group = []

# ----------------------------------- System Configuration general values
boxs = 1.43
#NBmethod = 'CutoffNonPeriodic' 
NBmethod = '' 
NBcut = 20
constraintTolerance = 0.0001
constraints = 'HBonds'
#rigidWater = True
#ewaldErrorTolerance = 0.0005
#hydrogenMass = 1.5*amu

# ----------------------------------- Simple parsing of args
#print("Full command:", sys.argv[0], sys.argv[1:])
i = 0 
for tok in sys.argv[1:20]:
      m = re.search(r'(?<=-)\w+', tok)
      if ( m ):
           i += 2 
           if ('ste' in m.group()):                # steps for min
               steps = int(sys.argv[i])
           elif (m.group() == 'af2'):              # AF2 imitate
               AF2 = True
               i -= 1
           elif (m.group() == 'af2w'):             # AF2 imitate plus implicit water
               AFw = True
               i -= 1
           elif (m.group() == 'rtt'):              # RTT method (more freedom)
               RTT = True
               i -= 1
           elif (m.group() == 'cif'):              # write in PDBx (CIF)
               mmCIF = True
               i -= 1
           elif ('rep' in m.group()):              # report on minimiztion steps 
               Reporting = True
               i -= 1
           elif ('gr' in m.group()):               # group all files at the end into a single file
               Grouping = True
               Groupfname = sys.argv[i]
           elif ('tol' in m.group()):              # tolerance for min
               tolerance = float(sys.argv[i])
           elif ('NBm' in m.group() ):             # NBmethod -- nonbondedMethod
               NBmethod = sys.argv[i]
           elif ('NBc' in m.group() ):
               NBcut = float(sys.argv[i])          # NBcutoff -- nonbondedCutoff
           elif (m.group() == 'inD'):
               indirname = sys.argv[i]
           elif (m.group() == 'outD'):
               outdirname = sys.argv[i]
           elif ('res' in m.group()):              # restrain atoms (move a little)
               atomsrestrained = (sys.argv[i])
           elif ('con' in m.group()):              # constrain atoms (no move at all)
               atomsconstrained = (sys.argv[i])
           elif (m.group() == 'kf'):
               stiffness =  float(sys.argv[i])
           elif (m.group() == 'box'):
               boxs = float(sys.argv[i])
           elif (m.group() == 'Ecut'):
               Ecut = float(sys.argv[i])
           elif ('lis' in m.group()):              # list
               with open(sys.argv[i], mode="r") as lg:
                    for ln in lg.readlines():
                        files.append(ln.strip().split()[0])
                        lEp.append(ln.strip().split()[0])
           elif ('for' in m.group() or 'ff' in m.group()):   # forcefield
               if ( sys.argv[i] == 'amber14'):
                  ffield = 'amber14-all.xml'
                  Water = False
               else:
                  ffield = 'amber99sb.xml'
                  Water = False
           elif ('wat' in m.group() ):              # water 
               if ( sys.argv[i] == 'amber14'):
                  Water = True
                  ffwater = 'amber14/tip3pfb.xml'   # tip3 real water
               elif ( sys.argv[i] == 'amber99'):
                  Water = True
                  ffwater = 'tip3p.xml'
               elif ( 'imp' in sys.argv[i] ):       # implicit water
                  Water = False
                  ffwater = 'implicit/gbn2.xml'
           else:
               print("Option %s not recognized" % tok)
               sys.exit(0)


# ---- Group all models into a single multimodel file
def GroupModels(Groupfname):
    OldF = ""
    nm = 0
    it = 0
    Ff = Files_to_group 
    with open(Groupfname, "w") as OA:
        for f in Ff:
            with open(f, "r") as IF:
                if f != OldF:
                    nm += 1
                    OldF = f
                    OA.write(f"MODEL      {nm}\n")
                    for line in IF:
                        if line.startswith("ATOM") or line.startswith("REMARK") or line.startswith("TER"):
                           OA.write(line)
                    OA.write("ENDMDL\n")
        OA.write("END\n")

# ----------------------------------- Define a subclass of Minimization Reporter
class Reporter(MinimizationReporter):
    interval = int 
    energies = []            # array to record progress
    def report(self, iteration, x, grad, args):
        # -- print current system energy to screen 
        if ( steps == 0 ):
    	    interval = 50   # report interval
        else:
            interval = steps/10 

        if iteration % interval == 0:
            print("    %12d  %12.4e  %12.4e" % (iteration, args['system energy']*toKcalmol, args['restraint energy']))

        # -- save energy at each iteration to an array we can use later
        self.energies.append(args['system energy'])

        # -- The report method must return a bool specifying if minimization should be stopped. 
        # -- You can use this functionality for early termination.
        return False

# -- Create an instance of our reporter
reporter = Reporter()

# ---------------------------------- if AF2 true mimetize all defaults as AF2 relax procedure
if AF2:
    Method = 'AF2' 
    Water = False
    ffield = 'amber99sb.xml'

    if steps == 99:            # adjust steps to 0 if not selected
       steps = 100
    else: 
       steps = steps 
    if ( tolerance == 99.0 ):  # adjust tolerance to 2.39 if not selected
       tolerance = 2.39 * ENERGY/(LENGTH)
    else:
       tolerance = tolerance * ENERGY/(LENGTH)
    if ( stiffness == 99.0 ):    # adjust stiffness, kf, if not selected 
       stiffness = 10.0 * ENERGY/(LENGTH**2)
    else:
       stiffness = stiffness * ENERGY/(LENGTH**2)
    if Ecut == 1.0e99 :            # adjust Ecutoff if not selected
       Ecut = 0.0
    else:
       Ecut = Ecut  
    if ( atomsrestrained == '' ):    # adjust restrained atoms
        atomsrestrained = ('heavy')
    if ( atomsconstrained != '' ):
        atomsrestrained = ''

    constraints = 'HBonds'
    integrator = openmm.LangevinIntegrator(0, 0.01, 0.0)
    NBmethod = 'NoCutoff' 
    NBcut = 999 
elif AFw:                   # -- AF2 plus implicit water, unlim steps, choice of tolerance, etc
    Method = 'AF2w' 
    Water = False
    ffield = 'amber99sb.xml'
    ffwater = 'implicit/gbn2.xml' 

    if steps == 99:            # adjust steps to 0 if not selected
       steps = 100
    else: 
       steps = steps 
    if ( tolerance == 99.0 ):  # adjust tolerance to 2.39 if not selected
       tolerance = 2.39 * ENERGY/(LENGTH)
    else:
       tolerance = tolerance * ENERGY/(LENGTH)
    if ( stiffness == 99.0 ):    # adjust stiffness, kf, if not selected 
       stiffness = 10.0 * ENERGY/(LENGTH**2)
    else:
       stiffness = stiffness * ENERGY/(LENGTH**2)
    if Ecut == 1.0e99 :            # adjust Ecutoff if not selected
       Ecut = 0.0
    else:
       Ecut = Ecut  
    if ( atomsrestrained == '' ):    # adjust restrained atoms
        atomsrestrained = ('heavy')
    if ( atomsconstrained != '' ):
        atomsrestrained = ''

    constraints = 'HBonds'
    integrator = openmm.LangevinIntegrator(0, 0.01, 0.0)
    NBmethod = 'NoCutoff' 
    NBcut = 999 
elif RTT:                   # -- RTT method 
    Method = 'RTT' 
    Water = False

    if steps == 99:            # adjust steps to 0 if not selected
       steps = 0
    else: 
       steps = steps 
    if ( tolerance == 99.0 ):  # adjust tolerance to 2.39 if not selected
       tolerance = 1.0 * ENERGY/(LENGTH)
    else:
       tolerance = tolerance * ENERGY/(LENGTH)
    if ( stiffness == 99.0 ):    # adjust stiffness, kf, if not selected 
       stiffness = 5.0 * ENERGY/(LENGTH**2)
    else:
       stiffness = stiffness * ENERGY/(LENGTH**2)
    if Ecut == 1.0e99 :            # adjust Ecutoff if not selected
       Ecut = 0.0
    else:
       Ecut = Ecut  
    if ( atomsrestrained == '' ):    # adjust restrained atoms
        atomsrestrained = 'backbone'
    if ( atomsconstrained != '' ):
        atomsrestrained = ''
    if ( ffield == '' ):             # adjust force field 
       ffield = 'amber99sb.xml'

    constraints = 'HBonds'
    integrator = openmm.LangevinIntegrator(0, 0.01, 0.0)
    NBmethod = 'NoCutoff' 
    NBcut = 99 
else:                       # -- everything else, we call it Custom
    Method = 'Custom' 
    Water = False    

    if steps == 99:            # adjust steps to 0 if not selected
       steps = 0
    else: 
       steps = steps 
    if ( tolerance == 99.0 ):  # adjust tolerance to 2.39 if not selected
       tolerance = 1.0 * ENERGY/(LENGTH)
    else:
       tolerance = tolerance * ENERGY/(LENGTH)
    if ( stiffness == 99.0 ):    # adjust stiffness, kf, if not selected 
       stiffness = 5.0 * ENERGY/(LENGTH**2)
    else:
       stiffness = stiffness * ENERGY/(LENGTH**2)
    if Ecut == 1.0e99 :            # adjust Ecutoff if not selected
       Ecut = 0.0
    else:
       Ecut = Ecut  
    # No atoms restrained nor constrained by default
    if ( ffield == '' ):             # adjust force field 
       ffield = 'amber99sb.xml'
    if ( NBmethod == '' ): 
       NBmethod = 'NoCutoff'
    NBcut = 99 

# ------------------------------------ final selection of nonbonded based of NBmethod
if ( NBmethod == 'NoCutoff' ):
       nonbondedMethod = NoCutoff 
elif ( NBmethod == 'CutoffPeriodic' ):
       nonbondedMethod = CutoffPeriodic 
elif ( NBmethod == 'CutoffNonPeriodic'):
       nonbondedMethod = CutoffNonPeriodic
elif ( NBmethod == 'PME' ):
       nonbondedMethod = PME 
elif ( NBmethod == 'Ewald' ):
       nonbondedMethod = Ewald 
else:
       print(" **err* NBmethod:", NBmethod, " not available")
       sys.exit(-1)

nonbondedCutoff = NBcut * LENGTH

# ------------------------------------ check whether directory already exists
# -- for outdir
if ( outdirname != './' ):
    currdir = os.getcwd()
    if not os.path.exists(currdir + '/' + outdirname):
       os.mkdir(currdir + '/' + outdirname)
       print("Folder %s created!" % (currdir + '/' + outdirname) )
    else:
       print("Folder %s already exists" % (currdir + '/' + outdirname))
       print("Pick another name, please")
       sys.exit(0)

# --- some header summarizing options being used
#print("\n ** OpenMM minimize: %d steps, tolerance %5.2f" % (steps,tolerance._value) )
print("\n **", sys.argv) 
print(" ** OpenMM minimize: ") 
print("             Method: %s" % Method )  
if steps > 0 :
   print("              Steps: %d" % steps )  
else: 
   print("              Steps: %d" % steps,"(Unlimited)" )  
print("          Tolerance: %.2f" % tolerance._value) 
print("         ForceField: %-20s %-10s" % (ffield, ffwater)) 

if ( atomsrestrained != '' ):
     print("   Atoms restrained: %-20s   Kf(stiffness): %s" % (atomsrestrained, stiffness)) 

if ( atomsconstrained != '' ):
     print("  Atoms constrained: %-20s" % (atomsconstrained)) 

if ( NBmethod == 'NoCutoff' ) :
     print("          Nonbonded: %-20s" % nonbondedMethod)
else:
     print("          Nonbonded: %-20s" % nonbondedMethod, "  Cut-off:", nonbondedCutoff)

print("           E cutoff: %6.1e \n" % Ecut) 

# ------------------------------------ cleaning of models with pdbfixer and then minimize it
if ( files == []):
       files = sys.argv[i+1:]

startTime0 = float(time.time())
for f in files:
       f = indirname + '/' + f 
       startTime = float(time.time())
       # ---------- preparing cosmetic file name removing dir ----------
       try: 
              idx = f.rindex('/')
              myp = f[0:idx+1]
       except: 
              myp = "" 
       name = f.removeprefix(myp)
       name = name.removesuffix('.pdb')
       name = name.removesuffix('_new')        # remove final _new if present 
       name = name.removeprefix('relaxed_')
       name = name.removeprefix('unrelaxed_') 

       # ------ fix structure with PDBfixer, add atoms as needed --------
       fixer = PDBFixer(filename=f)
       fixer.findMissingResidues()
       fixer.findNonstandardResidues()
       fixer.replaceNonstandardResidues()
       fixer.removeHeterogens(False)        # removes all HETATMs 
       fixer.findMissingAtoms()             # Finds ay missing heavy atom
       fixer.addMissingAtoms()              # will add heavy atoms for the missing residues
       fixer.addMissingHydrogens(7.0)       # will add hydrogens. pH = 7.0

       # ----------- add cubic box of water if requested -----------------
       if ( Water ):
            maxSize = max(max((pos[i] for pos in fixer.positions))-
                      min((pos[i] for pos in fixer.positions)) for i in range(3))
            boxSize = maxSize*Vec3(boxs, boxs, boxs)
            fixer.addSolvent(boxSize)

       # --------------------------------------------------- from PDB or mmCIF
       if mmCIF:
           fixed = outdirname+"/fixed_"+name+'.cif'        # output file for fixed coords
           outfile = outdirname+"/min_"+name+'.cif'	   # output file for minimized coords
       else:
           if AF2: 
              fixed = outdirname+"/fixed_"+name+'_af2.pdb' # output file for fixed coords
              outfile = outdirname+"/min_"+name+'_af2.pdb' # output file for minimized coords
           elif RTT:
              fixed = outdirname+"/fixed_"+name+'_rtt.pdb' # output file for fixed coords
              outfile = outdirname+"/min_"+name+'_rtt.pdb' # output file for minimized coords
           else:
              fixed = outdirname+"/fixed_"+name+'.pdb'     # output file for fixed coords
              outfile = outdirname+"/min_"+name+'.pdb'	   # output file for minimized coords
       # ---------------------------------------------------
       PDBFile.writeFile(fixer.topology, fixer.positions, open(fixed, 'w'))
       pdb = PDBFile(fixed)

       # ---------- adjust force field for explcit water if requested ----
       if ( Water ) : 
            forcefield = ForceField(ffield,ffwater)
       else:
            forcefield = ForceField(ffield)

       # ---------- Creation of the system based on coords from fixed PDB
       modeller = Modeller(pdb.topology, pdb.positions)
       #modeller.addHydrogens(forcefield,pH=7.0)
       system = forcefield.createSystem(modeller.topology, 
                                        nonbondedMethod=nonbondedMethod,
                                        nonbondedCutoff=nonbondedCutoff, 
                                        constraints=constraints)
       
       print("   ----   File:  %s \t (out: %s)" % (name, outfile)) 

       totalparticles = system.getNumParticles()

       # ----------- Restrain atoms in place if requested --------------------
       restrainedparticles = 0  
       if ( atomsrestrained != '' ):
            # -- for periodic boundadry conditions
            force = CustomExternalForce("0.5 * k * ((x-x0)^2+(y-y0)^2+(z-z0)^2)")
            force.addGlobalParameter("k", stiffness)
            system.addForce(force)
            for p in ["x0", "y0", "z0"]:
                force.addPerParticleParameter(p)
            for atom in pdb.topology.atoms():
                if atomsrestrained == 'heavy':
                    if (atom.element.name != 'hydrogen'):
                        force.addParticle(atom.index, pdb.positions[atom.index])
                        restrainedparticles += 1 
                elif atomsrestrained in ('backbone','bb') :
                    if (atom.name == 'N') or (atom.name == 'CA') or (atom.name == 'C'):
                        force.addParticle(atom.index, pdb.positions[atom.index])
                        restrainedparticles += 1 
                else:
                    if atom.name.casefold() in (atomsrestrained.casefold()):
                        force.addParticle(atom.index, pdb.positions[atom.index])
                        restrainedparticles += 1 

            print("    Restrained: %s (%d/%d particles)    k: %.2f %s" % 
                        (atomsrestrained,restrainedparticles,totalparticles,stiffness._value, stiffness.unit))

       # ---------- Constrain atoms to original position (no move) ------------
       constrainedparticles = 0 
       if ( atomsconstrained != '' ):
            for atom in pdb.topology.atoms():
                if atomsconstrained == 'heavy':
                    if (atom.element.name != 'hydrogen'):
                        system.setParticleMass(atom.index, 0)
                        constrainedparticles += 1
                elif atomsconstrained in ('backbone','bb') :
                    if (atom.name == 'N') or (atom.name == 'CA') or (atom.name == 'C'):
                        system.setParticleMass(atom.index, 0)
                        constrainedparticles += 1
                else:
                    if atom.name.casefold() in (atomsconstrained.casefold()):
                        system.setParticleMass(atom.index, 0)
                        constrainedparticles += 1
            print("      Constrained: %s (%d/%d particles) " %
                        (atomsconstrained,constrainedparticles, totalparticles))

       # for later getting the forces
       for i, frc in enumerate(system.getForces()):
           frc.setForceGroup(i)
       
       # -------- Build simulation context and minimize  --------------------
       if ( AF2 ):    # make sure if AF2 mimic use the LangevinIntegrator
           integrator = openmm.LangevinIntegrator(0, 0.01, 0.0)
       else:
           integrator = openmm.LangevinIntegrator(0, 0.01, 0.0)

       # --- intial state
       simulation = Simulation(modeller.topology, system, integrator)
       simulation.context.setPositions(modeller.positions)
       state = simulation.context.getState(getEnergy=True, getPositions=True)
       Eti  = state.getPotentialEnergy()
       posini = state.getPositions(asNumpy=True).value_in_unit(LENGTH)

       # ---- Minimize and report steps/Ep if requested
       if Reporting:
           print("            Iter    Epotential    Erestraint")
           print("          ------    ----------    ----------")
           simulation.minimizeEnergy(maxIterations=steps, tolerance=tolerance, reporter=reporter)
           print("          ----------------------------------")
       else:
           simulation.minimizeEnergy(maxIterations=steps, tolerance=tolerance)

       # --- final (minimized) state
       state = simulation.context.getState(getEnergy=True, getPositions=True)
       Etf  = state.getPotentialEnergy()
       positions = state.getPositions()
       posmin = state.getPositions(asNumpy=True).value_in_unit(LENGTH)

       # ---- evaluate RMSD for heavy 
       iat = 0
       itot = 0 
       diffsqr = 0.0 ;  
       for atom in pdb.topology.atoms(): 
           if (atom.element.name == 'hydrogen'):     # avoid H's always
               iat += 1 
           else: 
              for j in range(3): 
                 diffsqr += (posini[iat][j] - posmin[iat][j])**2
              itot += 1 
              iat += 1 
       rmsd = math.sqrt(diffsqr/itot)  
       # ---- as numpy thing, taking all atoms
       rmsd_all = np.sqrt(np.sum((posini - posmin)**2) / posini.shape[0])

       # ---- write down minimized file with Epot stored in REMARK record if
       #      the selected energy (Etf or Enb) is lower than cutoff (Ecut)

       # ---- Getting Enb (nombonded term) 
       for i, frc in enumerate(system.getForces()):
           state = simulation.context.getState(getEnergy=True, groups={i})
           if frc.getName() in 'NonbondedForce' : 
              Enb = state.getPotentialEnergy()._value 

       # -- write CIF and/or PDB formats if model gets accepted 
       Epot = f'{Etf._value:10.2f} J/mol'
       with open(outfile, mode="w") as file:
           if mmCIF : 
              if Enb < Ecut : 
                 PDBxFile.writeFile(simulation.topology, positions, file)
              else:
                 os.remove(outfile)
           else :
              if Enb < Ecut : 
                 if AF2: 
                    file.write("REMARK File: " + "min_"+name+"_af2" + "  Epot: " + Epot + "\n") 
                 elif RTT: 
                    file.write("REMARK File: " + "min_"+name+"_rtt" + "  Epot: " + Epot + "\n") 
                 else: 
                    file.write("REMARK File: " + "min_"+name+ "  Epot: " + Epot + "\n") 

                 PDBFile.writeFile(simulation.topology, positions, file)
                 Files_to_group.append(outfile)
              else:
                 os.remove(outfile)
              
       # -------- output summary info ------------------
       print("          Epot:   %12.4e kcal/mole   elapsed t: %12.4e s " % 
                (Etf._value*toKcalmol, float(time.time())-startTime))

       # get report on current forces
       print("\n          -- Sumary of forces")
       for i, frc in enumerate(system.getForces()):
           state = simulation.context.getState(getEnergy=True, groups={i})
           print("         %s \t %12.4e kcal/mol" % (frc.getName(), state.getPotentialEnergy()._value*toKcalmol))
       print("         rmsd: %.3f (%d atoms)" % (rmsd,itot) )
       print("         rmsd: %.3f (%d atoms)" % (rmsd_all,posini.shape[0]) )
       print(" ")
       # -------- remove intermediate fix PDB files (unrelaxed) ---
       os.remove(fixed)

# -------- group good mini models into a single file if requested
if Grouping:
   GroupModels(Groupfname)   

print("\n   Total time: ", str(datetime.timedelta(seconds=round(time.time()-startTime0))))
print('Done')

