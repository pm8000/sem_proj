import numpy as np
import matplotlib.pyplot as plt
import LB_trial_fun as fun
import LB_source as src

L=10 #domain length in lattice points
l=1 #pipe length
dx=l/L
dt=1e-6 #time step
t_end=8*dt #end time
nt=int(t_end/dt) #number of time steps
u_in=10 #inlet velocity
T_in=200+273 #inlet temperature
p_in=1e5 #inlet pressure
u_init=0 #initial velocity
p_init=1e5 #initial pressure
T_init=20+273 #initial temperature
cv=718 #specific heat coefficient for constant density (air)
gamma=1.4 #isentropic coefficient (air)
Pr=0.7 #Prandtl number
mu=1.5e-5 #dynamic viscosity (air)
d=0.01 #pipe diameter
l_f=0.026 #air thermal conductivity

R=cv*(gamma-1)
rho_in=p_in/(R*T_in) #inlet density
rho_init=p_init/(R*T_init) #initial density
E_init=cv*T_init + 0.5*u_init**2
alpha=src.get_alpha(u_in, d, mu/rho_in, Pr, l_f)

#non-dimensionalise
#C is used for normalising constants
#index nd means non dimensionalised

C_l=l
C_u=50*u_in
C_rho=rho_init
C_R=cv
C_T=C_u**2/C_R # not independant
cv_nd=cv/C_R
R_nd=R/C_R
mu_nd=mu/(C_rho*C_u*C_l)
dx_nd=dx/C_l
dt_nd=dt*C_u/C_l
t_end_nd=t_end*C_u/C_l
alpha_nd=alpha/(C_R*C_rho*C_u)


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
rho[:]=rho_init/C_rho
u[:]=u_init/C_u
T[:]=T_init/C_T
T_wall[:]=T_init/C_T
U[:]=(T_init/C_T*cv_nd)
E[:]=U[:]+u[:]**2/2

omega=fun.get_omega(T, R_nd, mu_nd, dt_nd)
fun.get_f_eq(f_eq, rho, u, T, R_nd, omega, dt_nd, dx_nd)
fun.get_g_eq(g_eq, rho, u, E, T, R_nd)
f[:,:]=f_eq[:,:]
g[:,:]=g_eq[:,:]

for t in range(nt):
    if t%1==0:
        print(t)
    #advection step
    f[0,:]=np.roll(f[0,:],-1)
    f[2,:]=np.roll(f[2,:], 1)
    g[0,:]=np.roll(g[0,:],-1)
    g[2,:]=np.roll(g[2,:], 1)

    #boundary condition
    f[0,-1]=0.001           #very small right propagating flow at the outlet
    f[2, 0]=u_in*rho_in/(C_u*C_rho)+f[0,0] #fix inlet density
    g[0,-1]=0.001           #very small right propagating flow at the outlet
    g[2, 0]=(cv_nd*T_in/C_T+0.5*(u_in/C_u)**2) - g[1,0] -g[0,0] #fix inlet energy

    #collision step
    rho=fun.get_rho(f)
    u=fun.get_u(f, rho)
    E=fun.get_E(g, rho)
    U=fun.get_U(E, u)
    T=fun.get_T(U, cv_nd)
    omega=fun.get_omega(T, R_nd, mu_nd, dt_nd)
    omega_1=fun.get_omega_1(omega, Pr)
    fun.get_f_eq(f_eq, rho, u, T, R_nd, omega, dt_nd, dx_nd)
    fun.get_g_eq(g_eq, rho, u, E, T, R_nd)
    fun.get_g_star(g_star, g_eq, g, rho, u, T, R_nd, dt_nd, dx_nd)
    src.get_source_populations(source, source_last, T, T_wall, alpha_nd, d/C_l, dt_nd, order=1)
    
    f[:,:]=omega*f_eq[:,:] + (1-omega)*f[:,:]
    g[:,:]=omega*g_eq[:,:] + (1-omega_1)*g[:,:] + (omega_1-omega)*g_star[:,:] + source[:,:]
    plt.plot(rho[:]*E[:], label=t)
plt.legend()
plt.show()
