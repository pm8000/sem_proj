
A clean version of the solver to get started with.
Modification in print statements and similar shall be commited to the rs_trial branch.
Run simulations by executing rs_trial.py, results are saved in the folder result (folder must exist before running)
For multithreading type mpirun -np num python rs_trial.py into your terminal and replace num by the amount of trheads
choose your parameters by editing the rs_trial.py file:
L=pipe length
dx=grid spacing
d=pipe diameter
t_end= simulation end time
cfl=cfl number
state= real or ideal (choose equation of state)
start_region= type 1 if initial condition is single phase, or 2 if initial condition is two phase, only considered for real gases
fluid=fluid name must be known for CoolProp (real gas only)
gamma=specific heat ratio, required for ideal gas only
R=specific ideal gas constant, required for ideal gas only
T_inlet=inlet temperature (K)
T_amb=initial temperature(K), considered if start_region=1
undercool= if start_region=2, the wall needs to be colder than boiling temperature to extract latent heat, specify here how much colder (K) than boiling temperature
q_init=initial steam quality, considered if start_region==2
p_amb=initial pressure (Pa)
u_inlet=inlet (and initial) velocity (m/s)
p_inlet=inlet pressure (Pa)
p_outlet=outlet pressure (Pa)

divisor= every n-th step will be written to the output file (if more than ca. 1,000,000 steps are written, the file cannot be opened by excell or similar)
appendix= part of the filename to specify the simulation
