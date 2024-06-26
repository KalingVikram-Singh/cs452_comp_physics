import numpy as np
import matplotlib.pyplot as plt
import math
from mpl_toolkits.mplot3d import Axes3D
def Ca(T):  #Specific heat capacity for copper lalttice
    # if T < 1358/20:
        return 0.36 + 8.6*10**(-5)*T + 2.9*10**(-9)*T**(2)
    # if T>1358/20:
    #     return 0.5
    
def Ce(T):  #Specific heat capacity for electronic lalttice
    # if T < 1795.5/20: 
        return 12.682 * T
    # if T > 1795.5/20 or T == 1795.5/20:
    #     return 227554.754
         
def Ka(T):   #Thermal heat conductuvity.
    # if T< 1358/20:
        return 0.6 + 0.0011*T - 2.6*10**(-7)*T**(2)
    # if T > 1358/20 :
        # return 2.1 
    
def dKadt(T):   #Derivative of thermal heat conductuvity with respect to temperature.
    # if T< 1358/20:
        return 0.0011 - 5.2*10**(-7)*T
    # if T > 1358/20 :
        # return 0.1

def dKedt(Tl):   #Thermal heat conductuvity for electronic lattifce.
    return -7 * 10**(2)/Tl

def Ke(Tl):   #Thermal heat conductuvity for electronic lattifce.
    return 3.5 * 10**(2)/(Tl**(2))

def A(i,k,dr,dt,n,Te):                                              # A term in the equation  -  Acts as source term for electrons, stopping power of electronic lattice
    alpha1 = 100
    alpha2 = 100
    t0 = 2.5*10**(-15) 
    r0 = 2.5*10**(-9)
    Se = 0.00000052
    T0 = 20
    Be = 1/((2*math.pi**(1.5))*(.84134*(r0**(2))*t0))
    A0 = Be*Se*alpha1*t0/(T0)
    resul = A0*(np.exp(-alpha2 * (i*dr)-(alpha1**(2) * (n*dt) - 0.05)**2 / 2))
    # print("A",resul)
    return resul

def B(i,k,dr,dt,n,Ta):                                              # B term in the equation  -  Acts as source term for ions, stopping power of atomic lattice
    alpha1 = 100
    alpha2 = 100
    r0 = 2.5*10**(-9)
    t0 = 2.5*10**(-15)                                            
    tau = 2.5*10**(-13) 
    R0 = 22*10**(-10)
    Sn = 0.00000086
    T0 = 20                                                        # Initial temperature, wihtout normalization
    Bn = 9*10**(-11)/((2*math.pi**(1.5))*(.84134*(R0**(2))*tau))
    B0 = Sn*Bn*alpha1*t0/(alpha2*r0*T0)
    # print(B0)
    # print("B = ",(B0 /(i*dr))*np.exp(-alpha2*r0*(dr*i) -n*dt*(alpha1*t0)/tau ))
    return B0*np.exp(-alpha2*r0*(dr*i) -n*dt*(alpha1*t0)/tau ) /(i*dr)


def heat_equation_explicit_integration( r, z, t, dt, dr, dz):                       # Function to solve the heat equation using explicit integration

    Te_initial = 1                                                                  # Initial temperature for electrons (normalized by room temperature (sample temperature))
    Ta_initial = 1                                                                  # Initial temperature for electrons (normalized by room temperature (sample temperature))
    g = 2.3 * 10**(10)                                                              # Electron-phonon coupling constant
    Te = np.zeros((r,z,t))
    Ta = np.zeros((r,z,t))                                                         # Creating grid for temperature loop
    alpha1 = 100                                                                   #Constants of transformation
    alpha2 = 100                                                                    #Constants of transformation
    t0 = 2.5*10**(-14)                                                              # Mean time for perpendicular electrons
    r0 = 2.5*10**(-9)                                                               # Space distance for highly excited electrons + ions

    # Set initial conditions
    Te[:,:,0] = Te_initial                                                          #At time = 0, the temperature is the initial temperature
    Ta[:,:,0] = Ta_initial
    Te[r-1,:,:] = Te_initial                                                        #At the boundaries, the temperature is the initial temperature
    Te[:,z-1,:] = Te_initial
    Ta[r-1,:,:] = Ta_initial
    Ta[:,z-1,:] = Ta_initial

    # Time stepping loop             
    for n in range(0,t-1):                      # This changes the time for in each step

        print("Time step: n = ",n,"\n")

        for i in range(1, r-1):                 # This is for the the radial distance
            # print("Spatial (r) step: i = ",i,"\n")

            for j in range(1, z-1):               # This is for the the axial distance

                # print("Spatial step (z): j = ",j,"\n")

                # Calculate spatial derivatives along with the initial conditions:
                if i==1:
                    dTe_dr = 0
                else:
                    dTe_dr = ((Te[i+1,j,n] - Te[i-1,j,n]) / (2*dr))                       # i -  R spatial variable

                if j==1:
                    dTe_dz = 0
                else:
                    dTe_dz = ((Te[i,j+1,n] - Te[i,j-1,n]) / (2*dz))                       # j -  Z spatial variable
                if i==1:
                    dTa_dr = 0
                else:
                    dTa_dr = ((Ta[i+1,j,n] - Ta[i-1,j,n]) / (2*dr))                       # i -  R spatial variable

                if j==1:
                    dTa_dz = 0
                else:
                    dTa_dz = ((Ta[i,j+1,n] - Ta[i,j-1,n]) / (2*dz))                     # j -  Z spatial variable

                d2Te_dr2 = ((Te[i+1,j,n] + Te[i-1,j,n] - 2*Te[i,j,n]) / (dr**2))
                d2Te_dz2 = ((Te[i,j+1,n] + Te[i,j-1,n] - 2*Te[i,j,n]) / (dz**2))
                d2Ta_dr2 = ((Ta[i+1,j,n] + Ta[i-1,j,n] - 2*Ta[i,j,n]) / (dr**2))
                d2Ta_dz2 = ((Ta[i,j+1,n] + Ta[i,j-1,n] - 2*Ta[i,j,n]) / (dz**2))


                # Update temperatures
                # print(d2Ta_dz2,d2Te_dz2,d2Ta_dr2,d2Te_dr2,dTe_dr,dTa_dr)
                # print((Ke(Ta[i,j,n])*dt)  *  (dTe_dr/(i*dr) + d2Te_dr2 + d2Te_dz2)/Ce(Te[i,j,n]) -  (dt/Ce(Te[i,j,n]))* (g*alpha1*t0)*(Te[i,j,n] - Ta[i,j,n]) + (dt/Ce(Te[i,j,n]))*A(i,j,dr,dt,n,Te),n)
                Te[i,j,n+1] = Te[i,j,n] + (Ke(Ta[i,j,n])*dt)  *  (dTe_dr/(i*dr) + d2Te_dr2 + d2Te_dz2)/Ce(Te[i,j,n]) -  (dt/Ce(Te[i,j,n]))* (g*alpha1*t0)*(Te[i,j,n] - Ta[i,j,n]) + (dt/Ce(Te[i,j,n]))*A(i,j,dr,dt,n,Te)
                Ta[i,j,n+1] = Ta[i,j,n] + (Ka(Ta[i,j,n])*dt)  *  (dTa_dr/(i*dr) + d2Ta_dr2 + d2Ta_dz2)/Ca(Ta[i,j,n]) +  (dt/Ca(Ta[i,j,n]))* (g*alpha1*t0)*(Te[i,j,n] - Ta[i,j,n]) + (dt/Ca(Ta[i,j,n]))*B(i,j,dr,dt,n,Ta)
    return Te, Ta
alpha1 = 100
t0 = 2.5*10**(-14)
dr = 10**(-3)                   # finite difference for radial direction
dz = 10**(-3)                   # finite difference for z direction
dt = 10**(-9)                # finite difference for time
r = 25
z = 25
t = 1000

print("The divisions in R, Z and T are:",r,z,t,"\n")

Te,Ta = (heat_equation_explicit_integration(r, z, t, dt, dr, dz))
R = np.linspace(0,r*dr,r)
Z = np.linspace(0,z*dz,z)
time = np.linspace(0,dr*t,t)

# Create a meshgrid and plot the graph.
import copy
total_temp =  copy.deepcopy(Te + Ta)

# print(total_temp)
R, Z = np.meshgrid(R, Z)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(R, Z, (total_temp[:,:,t-1]*20), cmap='viridis')
ax.set_xlabel('Radial Distance (x 250 nm)')
ax.set_ylabel('Axial Distance (x 250 nm)')
ax.set_zlabel('Temperature (K)')
plt.show()

# print(total_temp)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(R, Z, (Te[:,:,t-1]*20), cmap='viridis')
ax.set_xlabel('Radial Distance (x 250 nm)')
ax.set_ylabel('Axial Distance (x 250 nm)')
ax.set_zlabel('Temperature (K)')
plt.show()

# print(total_temp)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(R, Z, (Ta[:,:,t-1]*20), cmap='viridis')
ax.set_xlabel('Radial Distance (x 250 nm)')
ax.set_ylabel('Axial Distance (x 250 nm)')
ax.set_zlabel('Temperature (K)')
plt.show()

# Assuming R=10 and Z=20
R_index = 2
Z_index = 2

# Extract the temperature values for the specific R and Z
temp_values = (Te[R_index, Z_index, :]*20)

# Create an array for time
time = np.arange(0, len(temp_values))

# Plot Time against the temperature values
plt.plot(time, temp_values)
plt.xlabel(f'Time of evolution ( x {10**(-21)} s)')
plt.ylabel('Temperature (K)')
plt.show()
