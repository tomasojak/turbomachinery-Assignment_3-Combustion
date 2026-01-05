import numpy as np 
D=0.025 #m
U_u=30 #m/s
l_0=D/2 #integral length scale
u_0=U_u/20 #integral eddy velocity 

s_l0=2.2188282306851606 #m/s (laminar flame speed)
nu_u=2.0737213151400035e-05 #m/s^2 (bulk viscosity)

epsilon=np.power(u_0,3)/l_0 #rate of dissipation scale 
u_eta=np.power(nu_u*epsilon,0.25) #Kolmogorov velocity 

Ka=np.power(u_eta/s_l0,2) #Karlovitz number 

Re= U_u*D/nu_u

s=s_l0*np.power(1+np.power(u_0/s_l0,2),0.5) #flame speed 

print("Re", Re, "Ka", Ka, "s", s)