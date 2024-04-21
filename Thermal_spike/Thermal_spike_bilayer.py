import numpy as np
import matplotlib.pyplot as plt
import math
import copy
from mpl_toolkits.mplot3d import Axes3D
def Ca(T):  #Specific heat capacity for copper lalttice
        return 0.36 + 8.6*10**(-5)*T + 2.9*10**(-9)*T**(2)

    
def Ce(T):  #Specific heat capacity for electronic lalttice
        return 12.682 * T

         
def Ka(T):   #Thermal heat conductuvity.
        return 0.6 + 0.0011*T - 2.6*10**(-7)*T**(2)

    
def dKadt(T):   #Derivative of thermal heat conductuvity with respect to temperature.
        return 0.0011 - 5.2*10**(-7)*T

def dKedt(Tl):   #Thermal heat conductuvity for electronic lattifce.
    return -7 * 10**(2)/Tl

def Ke(Tl):   #Thermal heat conductuvity for electronic lattifce.
    return 3.5 * 10**(2)/(Tl**(2))

def Ca2(T):  #Specific heat capacity for copper lalttice
    if T < 1455/20:
        return 0.39 + 0.00019*T - 3.3*10**(-8)*T**(2) + 3.8*10**(-11)*T**(3)
    if T>1455/20:
        return 0.42
    
def Ce2(T):  #Specific heat capacity for electronic lalttice
    if T < 1795.5/20: 
        return 12.682 * T
    if T > 1795.5/20 or T == 1795.5/20:
        return 227554.754
         
def Ka2(T):   #Thermal heat conductuvity.
    if T< 1358/20:
        return 3.4 - 0.013*T - 2.2*10**(-5)*T**(2) - 1.5*10**(-8)*T**(3) + 3.6*10**(-12)*T**(4)
    if T > 1358/20 :
        return 0.5 
    
def dKadt(T):   #Derivative of thermal heat conductuvity with respect to temperature.
    if T< 1358/20:
        return 0.0011 - 5.2*10**(-7)*T
    if T > 1358/20 :
        return 0.1

def dKedt(Tl):   #Thermal heat conductuvity for electronic lattifce.
    return -7 * 10**(2)/Tl

def Ke2(Tl):   #Thermal heat conductuvity for electronic lattifce.
    return 3.5 * 10**(2)/(Tl**(2))




def A1(i,k,dr,dt,n,Te):                                              # A term in the equation  -  Acts as source term for electrons, stopping power of electronic lattice
    alpha1 = 100
    alpha2 = 100
    t0 = 2.5*10**(-15) 
    r0 = 2.5*10**(-9)
    Se = 0.00000049
    T0 = 300
    Be = 1/((2*math.pi**(1.5))*(.84134*(r0**(2))*t0))
    A0 = Be*Se*alpha1*t0/(T0)
    resul = A0*(np.exp(-alpha2 * (i*dr)-(alpha1**(2) * (n*dt) - 0.75)**2 / 2))
    # print("A",resul)
    return resul

def B1(i,k,dr,dt,n,Ta):                                              # B term in the equation  -  Acts as source term for ions, stopping power of atomic lattice
    alpha1 = 100
    alpha2 = 100
    r0 = 2.5*10**(-9)
    t0 = 2.5*10**(-15)                                            
    tau = 2.5*10**(-13) 
    R0 = 22*10**(-10)
    Sn = 0.00000090
    T0 = 20                                                       # Initial temperature, wihtout normalization
    Bn = 9*10**(-11)/((2*math.pi**(1.5))*(.84134*(R0**(2))*tau))
    B0 = Sn*Bn*alpha1*t0/(alpha2*r0*T0)
    # print(B0)
    # print("B = ",(B0 /(i*dr))*np.exp(-alpha2*r0*(dr*i) -n*dt*(alpha1*t0)/tau ))
    return B0*np.exp(-alpha2*r0*(dr*i) -n*dt*(alpha1*t0)/tau ) /(i*dr)
def A2(i,k,dr,dt,n,Te):                                              # A term in the equation  -  Acts as source term for electrons, stopping power of electronic lattice
    alpha1 = 100
    alpha2 = 100
    t0 = 2.5*10**(-15) 
    r0 = 2.5*10**(-9)
    Se = 0.000000175
    T0 = 300
    Be = 1/((2*math.pi**(1.5))*(.84134*(r0**(2))*t0))
    A0 = Be*Se*alpha1*t0/(T0)
    resul = A0*(np.exp(-alpha2 * (i*dr)-(alpha1**(2) * (n*dt) - 0.75)**2 / 2))
    # print("A",resul)
    return resul

def B2(i,k,dr,dt,n,Ta):                                              # B term in the equation  -  Acts as source term for ions, stopping power of atomic lattice
    alpha1 = 100
    alpha2 = 100
    r0 = 2.5*10**(-9)
    t0 = 2.5*10**(-15)                                            
    tau = 2.5*10**(-13) 
    R0 = 22*10**(-10)
    Sn = 0.000000867
    T0 = 20                                                       # Initial temperature, wihtout normalization
    Bn = 10**(-10)/((2*math.pi**(1.5))*(.84134*(R0**(2))*tau))
    B0 = Sn*Bn*alpha1*t0/(alpha2*r0*T0)
    # print(B0)
    # print("B = ",(B0 /(i*dr))*np.exp(-alpha2*r0*(dr*i) -n*dt*(alpha1*t0)/tau ))
    return B0*np.exp(-alpha2*r0*(dr*i) -n*dt*(alpha1*t0)/tau ) /(i*dr)
# Ideal contact: Let the first layer be 25 units long and second be same, 25 units long (6.25 nm each). For simplicity now, we take ideal contact and the material is copper, we can just compare the two effects.

def heat_equation_explicit_integration2layer( r, z, t, dt, dr, dz):                       # Function to solve the heat equation using explicit integration

    Te_initial = 1                                                                  # Initial temperature for electrons (normalized by room temperature (sample temperature))
    Ta_initial = 1                                                                  # Initial temperature for electrons (normalized by room temperature (sample temperature))
    g = 2.3 * 10**(10)                                                              # Electron-phonon coupling constant
    Te = np.zeros((r,z,t))
    Ta = np.zeros((r, z,t))                                                         # Creating grid for temperature loop
    Te2 = np.zeros((r,z,t))
    Ta2 = np.zeros((r, z,t))  

    alpha1 = 100                                                                   #Constants of transformation
    alpha2 = 100                                                                    #Constants of transformation
    t0 = 2.5*10**(-14)                                                              # Mean time for perpendicular electrons
    r0 = 2.5*10**(-9)                                                               # Space distance for highly excited electrons + ions

    # Set initial conditions
    Te[:,:,0] = Te_initial                                                          #At time = 0, the temperature is the initial temperature
    Ta[:,:,0] = Ta_initial
    
    Te2[:,:,0] = Te_initial                                                          #At time = 0, the temperature is the initial temperature for the second layer too.
    Ta2[:,:,0] = Ta_initial
    Te2[r-1,:,:] = Te_initial                                                        #At the boundaries, the temperature is the initial temperature
    Te2[:,z-1,:] = Te_initial
    Ta2[r-1,:,:] = Ta_initial
    Ta2[:,z-1,:] = Ta_initial


    # Time stepping loop             
    for n in range(0,t-1):                      # This changes the time for in each step

        print("Time step: n = ",n,"\n")
        # print(Te[4,4,n-1])

        for i in range(1, r-1):                 # This is for the the radial distance
            # print("Spatial (r) step: i = ",i,"\n")

            for j in range(1, z-1):               # This is for the the axial distance
                Te2[:,0,:] = Te[:,z-2,:] 
                Ta2[:,0,:] = Ta[:,z-2,:] 
                

                # print("Spatial step (z): j = ",j,"\n")

                # Calculate spatial derivatives along with the initial conditions:
                if i==1:
                    dTe_dr = 0
                if j < z-1:
                    dTe_dr = ((Te[i+1,j,n] - Te[i-1,j,n]) / (2*dr))                       # i -  R spatial variable
                elif j == z-1:
                    dTe_dr = ((Te[i+1,j,n] - Te[i-1,j,n]) / (2*dr))
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

                
                if j == 0:
                    dTa_dr2 = dTa_dr*Ka(Ta[i,j,n])
                    dTe_dr2 = dTe_dr*Ka(Ta[i,j,n])
                else:
                    dTa_dr2 = ((Ta[i+1,j,n] - Ta[i-1,j,n]) / (2*dr))
                    dTe_dr2 = ((Te[i+1,j,n] - Te[i-1,j,n]) / (2*dr))
                


                d2Te_dr2 = ((Te[i+1,j,n] + Te[i-1,j,n] - 2*Te[i,j,n]) / (dr**2))
                d2Te_dz2 = ((Te[i,j+1,n] + Te[i,j-1,n] - 2*Te[i,j,n]) / (dz**2))
                d2Ta_dr2 = ((Ta[i+1,j,n] + Ta[i-1,j,n] - 2*Ta[i,j,n]) / (dr**2))
                d2Ta_dz2 = ((Ta[i,j+1,n] + Ta[i,j-1,n] - 2*Ta[i,j,n]) / (dz**2))

                d2Te_dr22 = ((Te2[i+1,j,n] + Te2[i-1,j,n] - 2*Te2[i,j,n]) / (dr**2))
                d2Te_dz22 = ((Te2[i,j+1,n] + Te2[i,j-1,n] - 2*Te2[i,j,n]) / (dz**2))
                d2Ta_dr22 = ((Ta2[i+1,j,n] + Ta2[i-1,j,n] - 2*Ta2[i,j,n]) / (dr**2))
                d2Ta_dz22 = ((Ta2[i,j+1,n] + Ta2[i,j-1,n] - 2*Ta2[i,j,n]) / (dz**2))


                # Update temperatures
                
                Te[i,j,n+1] = Te[i,j,n] + (Ke(Ta[i,j,n])*dt)  *  (dTe_dr/(i*dr) + d2Te_dr2 + d2Te_dz2) -  (dt/Ce(Te[i,j,n]))* (g*alpha1*t0)*(Te[i,j,n] - Ta[i,j,n]) + (dt/Ce(Te[i,j,n]))*A1(i,j,dr,dt,n,Te)
                Ta[i,j,n+1] = Ta[i,j,n] + (Ka(Ta[i,j,n])*dt)  *  (dTa_dr/(i*dr) + d2Ta_dr2 + d2Ta_dz2) +  (dt/Ca(Ta[i,j,n]))* (g*alpha1*t0)*(Te[i,j,n] - Ta[i,j,n]) + (dt/Ca(Ta[i,j,n]))*B1(i,j,dr,dt,n,Ta)

                if j ==z-2:
                    Te2[:,0,n] = Ta[:,j,n]
                    Ta2[:,0,n] = Ta[:,j,n]

                Te2[i,j,n+1] = Te2[i,j,n] + (Ke2(Ta2[i,j,n])*dt)  *  (dTe_dr2/(i*dr) + d2Te_dr22 + d2Te_dz22) -  (dt/Ce2(Te2[i,j,n]))* (g*alpha1*t0)*(Te2[i,j,n] - Ta2[i,j,n]) + (dt/Ce2(Te2[i,j,n]))*A2(i,j,dr,dt,n,Te2)
                Ta2[i,j,n+1] = Ta2[i,j,n] + (Ka2(Ta2[i,j,n])*dt)  *  (dTa_dr2/(i*dr) + d2Ta_dr22 + d2Ta_dz22) +  (dt/Ca2(Ta2[i,j,n]))* (g*alpha1*t0)*(Te2[i,j,n] - Ta2[i,j,n]) + (dt/Ca2(Ta2[i,j,n]))*B2(i,j,dr,dt,n,Ta2)
        Te[:,z-1,t-1] = 0.7*Te[:,z-2,t-1]
        Ta[:,z-1,t-1] = 0.7*Ta[:,z-2,t-1]
    return Te, Ta ,Te2, Ta2

alpha1 = 100
t0 = 2.5*10**(-14)
dr = 10**(-3)                   # Finite difference for radial direction.
dz = 10**(-3)                   # Finite difference for z direction.
dt = 10**(-9)                   # Finite difference for time.
r = 25
z = 25
t = 8500
print("The divisions in R, Z and T are:",r,z,t,"\n")

Te,Ta,Te2,Ta2 = heat_equation_explicit_integration2layer(r, z, t, dt, dr, dz)
Te_1 = copy.deepcopy(Te)
Ta_1 = copy.deepcopy(Ta)
Te_2 = copy.deepcopy(Te2)
Ta_2 = copy.deepcopy(Ta2)
R = np.linspace(0,r*dr,r)
Z = np.linspace(0,z*dz,z)
time = np.linspace(0,dr*t,t)

# Create a meshgrid and plot the graph.
total_temp =  copy.deepcopy(Te + Ta)
total2_temp = copy.deepcopy(Te2 + Ta2)
net1 = np.dstack((Te, Te2))                             
nat2 = np.dstack((Ta, Ta2))
system_temp = np.dstack((total_temp, total2_temp))

# Assuming R=16 and Z=20
R_index = 16
Z_index = 20
R, Z = np.meshgrid(R, Z)                                                        #Making a grid for the plotting in 3

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_title("Total temperature of system") 
ax.plot_surface(R, Z, (total_temp[:,:,t-1]*20), cmap='viridis')          #Gives total temperature vs space coordinates
ax.set_xlabel('Radial Distance (x 250 nm)')
ax.set_ylabel('Axial Distance (x 250 nm)')
ax.set_zlabel('Temperature (K)')
plt.show()

# print(total_temp)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_title("Lattice temperature of system") 
ax.plot_surface(R, Z, (Ta[:,:,t-1]*20), cmap='viridis')                  #Gives Lattice temperature vs space coordinates
ax.set_xlabel('Radial Distance (x 250 nm)')
ax.set_ylabel('Axial Distance (x 250 nm)')
ax.set_zlabel('Temperature (K)')
plt.show()


# Extract the temperature values for the specific R and Z
temp_values = (total_temp[R_index, Z_index, :]*20)

# Create an array for time
time = np.arange(0, len(temp_values))

# Plot Time against the temperature values
plt.plot(time, temp_values)
plt.xlabel(f'Time of evolution ( x {10**(-21)} s)')                     #Gives  temperature vs time
plt.ylabel('Temperature (K)')
plt.show()

# Create 2*R and 2*Z
R = np.linspace(r*dr, 2*r*dr, r)
Z = np.linspace(z*dz, 2*z*dz, z)
R, Z = np.meshgrid(R, Z)  

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_title("Lattice temperature of system 2") 
ax.plot_surface(R, Z, (Ta2[:,:,t-1]*20), cmap='viridis')                 #Gives elattice temperature vs space coordinates of 2nd system
ax.set_xlabel('Radial Distance (x 250 nm)')
ax.set_ylabel('Axial Distance (x 250 nm)')
ax.set_zlabel('Temperature (K)')
plt.show()



# Create 2*R and 2*Z
R = np.linspace(0, 2*r*dr,r)
Z = np.linspace(0, 2*z*dz,z)

# Create a meshgrid
R, Z = np.meshgrid(R, Z)                                               #Making a grid for the plotting in 3D for 2 system


nat2 = np.concatenate((Ta_1,Ta_2),axis = 2)
# Plot the surface
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_title("Lattice temperature of total system 1 and system 2") 
ax.plot_surface(R, Z, (nat2[:,:,t-1]*20), cmap='viridis')
ax.set_xlabel('Radial Distance (x 250 nm)')
ax.set_ylabel('Axial Distance (x 250 nm)')
ax.set_zlabel('Temperature (K)')
plt.show()


# Ensure that the third argument in np.linspace matches the shape of temp_values
R = np.linspace(0, 2*r*dr)
Z = np.linspace(0, 10*z*dz,2*z)
print(R)


print(Z)
nat2 = np.concatenate((Ta_1,Ta_2),axis = 1)
print(len(Z))
print(len(nat2[10, :, t-1]))
# Extract the temperature values for the specific R and Z
temp_values = (nat2[10, :, int((t-1))]*20)
# Create an array for time
# Plot Time against the temperature values
plt.plot(Z, temp_values,label = f"Temperature profile at r = 10 x2.5Angs")

temp_values = (nat2[5, :, int((t-1))]*20)
# Create an array for time
# Plot Time against the temperature values
plt.plot(Z, temp_values,label = f"Temperature profile at r = 5 x2.5Angs")

temp_values = (nat2[20, :, int((t-1))]*20)
# Create an array for time
# Plot Time against the temperature values
plt.plot(Z, temp_values,label = f"Temperature profile at r = 20 x2.5Angs")



plt.legend()
plt.xlabel(f' Axial distance ( x 50 Angs s)')                     #Gives  temperature vs time
plt.ylabel('Temperature (K)')
plt.xlim(0,0.233)
plt.show()
print(Te2[10,0,t-1])
