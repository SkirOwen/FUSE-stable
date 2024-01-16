from fuse202.all import *

run_fuse(

# FUSE setups ##################################################################
composition={'Ca':3,'Ti':2,'O':7}, # the emperical formula unit for this calculation.
max_atoms=50,  # Maximum number of atoms to use in the overall calculation
imax_atoms=50, # Maximum number of atoms to use in only the initial population, for large systems, it can be helpful to set this smaller than the "max_atoms".
restart=False,  # restart a previous calculation? If set to True, FUSE will attempt to restart from a previous calculation.
read_exisiting_structures = True, # read in and slice previously generated structures for modules? When using the ML structure geneartaion, set this to True.
path_to_structures = 'reference_structures', # should be a directory containing structures from the ML model, this should be left as this value. This directory should exist prior to starting the calculation. In a future version of the code we will get FUSE to create this directory automatically!
initial_gen=25, # the number of structures to include in the initial population
iterations=15,  # the number of structures to perform energy calculations on in this run of FUSE.

################################################################################

# Specific bits for the BH search rountine #####################################

################################################################################

search_gen_bh=1, # Number of new structures to generate at each step in the basin hopping search
melt_threshold=150, # After this number of steps since the bashin hopping routine makes a downhill step, FUSE will increase the temperateure parameter to escape local minima
rmax=500, # The basin hopping routine will be considered converged if this many structures are generated since the current lowest energy structure was located.
T=0.02, # the temperature parameter for the basin hopping routine, the higher the value, the larger uphill step that the basin hopping routine will accept.

################################################################################

# Specific bits for the Gn-Boss ML model for structure generation ##############

################################################################################

#variables used to run ML structure generation:
#gn boss model for structure generation
generate_gn_boss_structures=True, # if set to true, when FUSE is firt launched, it will run gn-boss to generate the pool of referennce structures for this calculation.
gn_boss_command= '$CONDA/conda run -p $GNBOSS/ python get_cifs_for_FUSE.py', # For my machine, I've setup gn-boss in a seperate python environment, this is the command for that version of python. You should not need to edit this if you have set the enviroment variables defined in the README file for FUSE v.2
gn_search='tpe', # 'rand' random search, 'tpe' baysian opt, 'pso' particle swarm
gn_max_step=5000, # number of generation attempts for gn-boss
gn_template_path= os.environ['GNBOSS_TEMP'], #path to template files for using gn-boss You should not need to edit this if you have set the enviroment variables defined in the README file for FUSE v.2 
gn_zn_range=[1,2,3,4], # numbers of formula units to scan with GN-BOSS for generating structures
rank_gn_structures='opti', #if None; do not rank structures, this should only be set if pull_random = True above, if 'opti' rank with SPPs AFTER geometry optimising the, if 'sing' rank based on single point calculations with SPPs. 
clear_previous_gn_structures=True, #if set to True, before starting the calcluation, remove any previous structures from reference structures & gn-boss generated_results.
generate_structures_only=False, # Run the ML structure generation & ranking, then exit if set to True. To then re-start with the basin hopping routine, set restart above to False along with generate_gn_boss_structures and clear_previous_gn_structures. 

# GULP options to rank the structures generated by the ML model above.
# For more details on how to use GULP with ASE please see the documentation here:
# https://wiki.fysik.dtu.dk/ase/ase/calculators/gulp.html

r_kwds=['opti conj conp noelectro','opti conj conp noelectro','opti conp noelectro','opti conp noelectro','sing conp noelectro'], # keywords for the gulp input, 
#multiple strings indicate running GULP in multiple stages
r_gulp_opts=[
['\ninclude ./lib.lib\ndump temp.res\nmaxcyc 250\nstepmax 0.001\ntime 5 minutes\ngmax 0.1\ngtol 0.1\nftol 0.1\nxtol 0.1'],
['\ninclude ./lib.lib\ndump temp.res\nmaxcyc 250\nstepmax 0.5\ntime 5 minutes\ngmax 0.1\ngtol 0.1\nftol 0.1\nxtol 0.1'],
['\ninclude ./lib.lib\ndump temp.res\nmaxcyc 1500\nlbfgs_order 5000\ntime 5 minutes\ngmax 0.1\ngtol 0.1\nftol 0.1\nxtol 0.1'],
['\ninclude ./lib.lib\ndump temp.res\nmaxcyc 1500\nlbfgs_order 5000\ntime 5 minutes\ngmax 0.1\ngtol 0.1\nftol 0.1\nxtol 0.1'],
['\ninclude ./lib.lib\ndump temp.res\nmaxcyc 1500\nlbfgs_order 5000\ntime 5 minutes\ngmax 0.1\ngtol 0.1\nftol 0.1\nxtol 0.1'],
],	# options for gulp, must include one line per set of inputs in "kwds"
r_lib='dummy.lib', # library file for interatomic potentials	

################################################################################

# Definitions for energy calculator(s) #########################################

################################################################################

### gulp inputs
ctype='mixed', # flag to pass to ase inorder to use gulp as the energy minimiser

calcs=['gulp','vasp'],
# Options for VASP

# For more information to configure Vasp, please see the ase documentation available here:
# https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#vasp
vasp_opts={
'1':Vasp(prec='Fast',ediff=1.E-4,nelm=50,encut=520,ibrion=2,isif=3,nsw=100,ediffg=0.05,nwrite=1,ncore=18,algo='Fast',gamma=True,lreal='Auto',setups={'Ca':'_sv'},kspacing=0.5),
'2':Vasp(prec='Fast',ediff=1.E-4,nelm=50,encut=520,ibrion=2,isif=3,nsw=50,ediffg=0.05,nwrite=1,ncore=18,algo='Fast',gamma=True,lreal='Auto',setups={'Ca':'_sv'},kspacing=0.5),
'3':Vasp(prec='Normal',ediff=1.E-5,nelm=50,encut=520,ibrion=1,isif=3,nsw=10,ediffg=-0.03,nwrite=1,ncore=18,gamma=True,lreal='.False.',setups={'Ca':'_sv'},kspacing=0.3,ismear=0),
'4':Vasp(prec='Normal',ediff=1.E-5,nelm=50,encut=520,ibrion=1,isif=3,nsw=10,ediffg=-0.03,nwrite=1,ncore=18,gamma=True,lreal='.False.',setups={'Ca':'_sv'},kspacing=0.3,ismear=0),
'5':Vasp(prec='Normal',ediff=1.E-5,nelm=50,encut=520,ibrion=1,isif=3,nsw=50,ediffg=-0.03,nwrite=1,ncore=18,gamma=True,lreal='.False.',setups={'Ca':'_sv'},kspacing=0.3,ismear=0),
},
kcut=[20],

#Options for GULP
# For more details on how to use GULP with ASE please see the documentation here:
# https://wiki.fysik.dtu.dk/ase/ase/calculators/gulp.html

kwds=['opti conp conj noelectro','opti conp noelectro conj','opti conp conj noelectro','opti conp conj noelectro','sing conp noelectro'], # keywords for the gulp input,

gulp_opts=[
['\ndump temp.res\nmaxcyc 250\ntime 5 minutes\ncutp 11.0 2.0\ninclude ./lib.lib\nstepmax 0.001'],
['\ndump temp.res\nmaxcyc 250\ntime 5 minutes\ncutp 11.0 2.0\ninclude ./lib.lib\ngmax 0.1\ngtol 0.1\nftol 0.1\nxtol 0.1\stepmax 0.5'],
['\ndump temp.res\nmaxcyc 1500\ntime 5 minutes\ncutp 11.0 2.0\ninclude ./lib.lib\ngmax 0.1\ngtol 0.1\nftol 0.1\nxtol 0.1'],
['\ndump temp.res\nmaxcyc 1500\ntime 5 minutes\ncutp 11.0 2.0\ninclude ./lib.lib\ngmax 0.1\ngtol 0.1\nftol 0.1\nxtol 0.1\nswitch_minimiser bfgs gnorm 1'],
['\ndump temp.res\nmaxcyc 1500\ntime 5 minutes\ncutp 11.0 2.0\ninclude ./lib.lib\ngmax 0.1\ngtol 0.1\nftol 0.1\nxtol 0.1\nswitch_minimiser bfgs gnorm 1'],

],   # options for gulp, must include one line per set of inputs in "kwds"
lib='dummy.lib', # library file for interatomic potentials, this is a dummy file, and needs to exist prior to starting the calculation. This will be fixed in a future version.
shel=[''],   # species using shells in gulp

# GULP options for FUSE:

gulp_timeout='', # timeout for gulp calls in seconds, leave as '' to ignore, use on Windows machines only.
assemble_spp_=True, # create a fresh SPP potential for this calculation, for most applications, this should be set to True.
spp_path=os.environ['SPP_PATH'], # if you have set teh enviroment variables required in the README file for FUSE v.2, you should not need to change this.
gulp_command=os.environ['ASE_GULP_COMMAND'], # if you have set teh enviroment variables required in the README file for FUSE v.2, you should not need to change this.


)
