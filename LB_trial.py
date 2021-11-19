import numpy as np
import matplotlib.pyplot as plt
import LB_trial_fun as fun
import LB_source as src

L=100 #domain length in lattice points
dt=1 #time step
dx=1 #spacial spacing
t_end=5e3 #end time
nt=int(t_end/dt) #number of time steps
u_in=10 #inlet velocity
T_in=200 #inlet temperature
p_in=1e5 #inlet pressure
u_init=0 #initial velocity
p_init=1e5 #initial pressure
T_init=20+273 #initial temperature
cv=718 #specific heat coefficient for constant density (air)
gamma=1.4 #isentropic coefficient (air)
Pr=0.7 #Prandtl number
mu=1.5e-5 #dynamic viscosity (air)
d=0.01 #pipe diameter
l_f=0.026 #air

R=cv*(gamma-1)
rho_in=p_in/(R*(T_in+273)) #inlet density
rho_init=p_init/(R*(T_init+273)) #initial density
E_init=cv*T_init + 0.5*u_init**2
alpha=src.get_alpha(u_in, d, mu/rho_in, Pr, l_f)

f=np.zeros([3,L]) #first row with velocity -c, second row with velocity 0 and third row with velocity c
f_eq=np.zeros([3,L])
g=np.zeros([3,L])
g_eq=np.zeros([3,L])
g_star=np.zeros([3,L])
source=np.zeros([3,L])
source_last=np.zeros(L)
rho=np.zeros(L) #density
u=np.zeros(L) #velocity
E=np.zeros(L) #total energy
U=np.zeros(L) #internal energy
T=np.zeros(L) #temperature
T_wall=np.zeros(L) #wall temperture
omega=np.zeros(L)
omega_1=np.zeros(L)
#initial condition
rho[:]=rho_init
u[:]=u_init
T[:]=T_init
T_wall[:]=T_init
U[:]=(T_init*cv)
E[:]=U[:]+u[:]**2/2

omega=fun.get_omega(T, R, mu, dt)
fun.get_f_eq(f_eq, rho, u, T, R, omega, dt, dx)
fun.get_g_eq(g_eq, rho, u, E, T, R)
f[:,:]=f_eq[:,:]
g[:,:]=g_eq[:,:]

for t in range(nt):
    #advection step
    f[0,:]=np.roll(f[0,:],-1)
    f[2,:]=np.roll(f[2,:], 1)
    g[0,:]=np.roll(g[0,:],-1)
    g[2,:]=np.roll(g[2,:], 1)

    #boundary condition #TODO
    f[0,-1]=0.1
    f[2, 0]=0.1
    g[0,-1]=0.1
    g[2, 0]=0.1

    #collision step
    rho=fun.get_rho(f)
    u=fun.get_u(f, rho)
    E=fun.get_E(g, rho)
    U=fun.get_U(E, u)
    T=fun.get_T(U, cv)
    omega=fun.get_omega(T, R, mu, dt)
    omega_1=fun.get_omega_1(omega, Pr)
    fun.get_f_eq(f_eq, rho, u, T, R, omega, dt, dx)
    fun.get_g_eq(g_eq, rho, u, E, T, R)
    fun.get_g_star(g_star, g_eq, g, rho, u, T, R, dt, dx)
    src.get_source_populations(source, source_last, T, T_wall, alpha, d, dt)
    
    f[:,:]=omega*f_eq[:,:] + (1-omega)*f[:,:]
    g[:,:]=omega*g_eq[:,:] + (1-omega_1)*g[:,:] + (omega_1-omega)*g_star[:,:] + source[:,:]
