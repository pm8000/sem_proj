import numpy as np

def get_rho(f):
    #calculate macroscopic local density
    #intput: population matrix (3xN)
    #output: density vector (1xN)
    return f[0,:]+f[1,:]+f[2,:]

def get_u(f, rho):
    #calculate macroscopic local velocity
    #intput: population matrix f (3xN), macroscopic local density rho
    #output: velocity vector (1xN)
    return (f[2,:]-f[0,:])/rho[:]

def get_E(f, rho):
    #calculate macroscopic local total energy
    #input: population matrix f (3xN), macroscopic local density rho
    #output: total energy vector (1xN)
    return (f[0,:]+f[1,:]+f[2,:])/(rho[:])

def get_U(E, u):
    #calculate internal energy
    #input:  (1xN) vectors of macroscopic local velocity (u) and total energy (E)
    #output: internal energy vector (1xN)
    return E[:]-0.5*u[:]**2

def get_T(U, cv):
    #calculate temperature field
    #intput: (1xN) vector of internal energy field
    #output: (1xN) temperature field vector
    return U[:]/cv

def get_P(f):
    #calculate total pressure
    #input: population matrix (3xN)
    #output: local total pressure field vector (1xN)
    return f[0,:]+f[2,:]

def get_P_eq(rho, u, T, R):
    #calculate equilibrium total pressure
    #input: vectors (1xN) of density (rho) and velocity(u) and temperatrue (T) field, scalar R (gas constant)
    #output: 1xN vector with equilibrium total pressure
    return rho[:]*(u[:]**2 + R*T)

def get_Q(rho, u, T, R):
    #get anisotropy correction for third order moment
    #input: vectors (1xN) of density (rho) and velocity(u) and temperatrue (T) field, scalar R (gas constant)
    #output: 1xN vector with anisotropy correction
    return rho[:]*u[:]*(1-3*R*T[:])-rho[:]*u[:]**3

def get_dQdx(Q, dx):
    #get first derrivative of Q
    #a second order central scheme is used, gradients are assumed to be zero otside of domain
    #input: (1xN) vector of Q and grid spacing dx
    #output: (1xN) vector of first derrivative of Q
    dQdx=(np.roll(Q, 1)[:]-np.roll(Q,-1)[:])/(2*dx)
    #apply BC
    dQdx[0]=(Q[1]-Q[0])/(2*dx)
    dQdx[-1]=(Q[-1]-Q[-2])/(2*dx)
    return dQdx

def get_P_ex(P_eq, dQdx, rho, dt, omega):
    #compensate equilibrium pressure
    #inputs: (1xN) vectors of equilibrium pressure (P_eq), derrivative of Q (dQdx), density (rho) and scalars omega and dt
    #output: (1xN) vector of compensated pressure
    return P_eq[:] + dt*(2-omega)/(2*rho[:]*omega)*dQdx[:]

def get_f_eq(f_eq, rho, u, T, R, omega, dx, dt):
    #calculate equilibrium populations for momentum lattice
    #input: solution matrix f_eq (3xN), vectors (1xN) for density (rho), velocity (u), temperature (T) and scalars R, omega, dx, dt
    #output: the equilibrium populations will be stored in f_eq (3xN matrix), passed by reference
    P_eq=get_P_eq(rho, u, T, R)
    Q=get_Q(rho, u, T, R)
    dQdx=get_dQdx(Q, dx)
    P_ex=get_P_ex(P_eq, dQdx, rho, dt, omega)
    f_eq[0,:]=0.5 * rho[:] * (-u[:] + P_ex[:])
    f_eq[1,:]=rho[:] * (1 - P_ex[:])
    f_eq[2,:]=0.5 * rho[:] * (u[:] + P_ex[:])

def get_g_eq(g_eq, rho, u, E, T, R):
    #calculate equilibrium population for energy
    #inputs: solution matrix g_eq (3xN), vector (1xN) for density (rho), velocity (u), total energy (E) and temperatrue (T) field, scalar R
    #output: the equilibrium populations will be stored in f_eq (3xN matrix), passed by reference
    P_eq=get_P_eq(rho, u, T, R)
    q_eq=(E[:] + R*T[:]) * rho[:] * u[:]
    R_eq=(E[:] + R*T[:]) * P_eq + rho[:]*R*T[:]*u[:]**2
    g_eq[0,:]=0.5 * (-q_eq[:] + R_eq[:])
    g_eq[1,:]=rho[:] * E[:] - R_eq
    g_eq[2,:]=0.5 * (q_eq[:] + R_eq[:])

def get_g_star(g_star, g_eq, g, rho, u, T, R, dt, dx):
    #calculate quasiequilibrium population for energy
    #inputs: solution matrix g_star (3xN), matrices (3xN) for equilibrium populations of energy (g) and energy populations (g)
    #           vectors (1xn) for density (rho), velocity (u) and temperatrue (T)
    #           scalars R, dt, dx
    #output: the quasiequilibrium populations will be stored in f_eq (3xN matrix), passed by reference
    P=get_P(g)
    P_eq=get_P_eq(rho, u, T, R)
    Q=get_Q(rho, u, T, R)
    dQdx=get_dQdx(Q, dx)
    g_star[0,:]=g_eq[0,:] - 0.5 * u[:] * (P[:] - P_eq[:] + dt/2*dQdx[:])
    g_star[1,:]=g_eq[1,:]
    g_star[2,:]=g_eq[2,:] + 0.5 * u[:] * (P[:] - P_eq[:] + dt/2*dQdx[:])

def get_omega(T, R, mu, dt):
    #calculate relaxation parameter
    #input: solution vector (1xN) omega, (1xN) temperature field T, scalars R, mu, dt
    #output: (1xN) vector omega
    return 2*R*dt*T[:]/(2*mu+R*T[:]*dt)

def get_omega_1(omega, Pr):
    #calculate relaxation parameter for quasiequilibrium
    #input: solution vector (1xN) omega_1, relaxation parameter (1xN) omega, scalar Prandtl mumber (Pr)
    #output: quasiequilibrium relaxation parameter
    return 2*omega[:]*Pr/(omega[:]*(Pr-1)+2)
