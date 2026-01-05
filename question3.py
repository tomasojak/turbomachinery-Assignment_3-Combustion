import numpy as np 
D=0.025 #m
l_0=D/2 #integral length scale

s_l0=2.2188282306851606 #m/s (laminar flame speed)
nu_u=2.0737213151400035e-05 #m/s^2 (bulk viscosity)

U_u_fb=np.power(400/19/21,0.5)*s_l0

print(U_u_fb, "U_u,fb")

Re=U_u_fb*D/nu_u
u_0=U_u_fb/20 #integral eddy velocity 
epsilon=np.power(u_0,3)/l_0 #rate of dissipation scale 
u_eta=np.power(nu_u*epsilon,0.25) #Kolmogorov velocity 

Ka=np.power(u_eta/s_l0,2) #Karlovitz number 

print("Ka", Ka, "Re", Re)