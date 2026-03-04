import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Parameters
Lambda_h = 5      # Human recruitment rate
mu_h = 1/(60*12)   # Human natural death rate (per month)
mu_v = 1.87   # Vector natural death rate (converted to per month)
nu_h = 0.05
omega_h=1/60       # Vaccination rate for humans
alpha_h = 1/3      # Immunity loss rate for humans (per month)
gamma_h = 30/14    # Progression rate from exposed to infected for humans (per month)
gamma_v = 30/14  # Progression rate for vectors (converted to per month)
p = 0.4            # Proportion progressing to infectious stage in humans
sigma_h = 1/24     # Recovery rate for asymptomatic humans (per month)
kappa_h = 0.1      # Rate of progression from asymptomatic to infectious
delta_h = 0.00274  # Disease-induced death rate for infectious humans
rho_vh = 0.22       # Transmission rate from vectors to humans
rho_hv = 0.48       # Transmission rate from humans to vectors
bite_rate = 15      # Maximal biting rate
K_b = 60         # Saturation constant for biting rate
B = 30
r=50               # Initial number of breeding sites
irrigation_factor = 1.2  # Irrigation multiplier affecting breeding sites
tau_max = 10/95    # Max recovery rate near health centers
lambda_base = 1  # Base infection rate
d_bs = 1           # Distance to breeding sites (affects susceptibility)
d_hc = 1          # Distance to health center (affects recovery)
f_s, c_s, f_a, c_a, f_i, c_i, b_e = 0, 0.5, 0, 0.5, 0, 0.5, 0.5

# Function to model seasonality using the sinusoidal function
def seasonality(t):
    return 1.0 + 0.5* np.sin((2 * np.pi / 12) * (t - 6))#+ 0.1* np.cos((2 * np.pi / 12) * (t - 7))
    #return 2.4172+ 1.4061*np. sin(np.pi*t/6) + 1.4061*np. cos(np.pi*t/6) - 0.0817*np. cos(np.pi*t/3) + 0.0969*np.cos(np.pi*t/2)-0.0294*np. cos(2*np.pi*t/3) + 0.0379*np.cos(5*np.pi*t/6) - 0.7357*np.sin(np.pi*t/6) - 0.0249*np.sin(np.pi*t/3)+0.1086*np.sin(np.pi*t/2) + 0.0012*np.sin(2*np.pi*t/3) - 0.001*np.sin(5*np.pi*t/6) 
# Biting rate function with breeding site and seasonal variation
def biting_rate(t):
    return bite_rate*seasonality(t)

def nu_h_t(t):
    if t <= 721:
        return 0  # No vaccination for the first 120 months
    else:
        return 0.05  # Vaccination starts after 120 months


# Bed net usage during rainy months (increases during months 7-10)
def bed_net_usage(t):
    month = t % 12
    if 7 <= month <= 10:
        return 0.6  # Increased efficacy during rainy months
    else:
        return 0.3  # Default efficacy
    
# Determine if chemotherapy is being applied during the current month
def P(t, A_h, I_h):
    month = t % 12
    if t >= 601 and 7 <= month <= 10:
        p = 0.25 * 0.5  # Reduce sigma during chemotherapy (after 120 months)
    else:
        p = 0.5  # Normal sigma without chemotherapy or before 120 months
    return p
    

# Recovery rate based on distance to health center
def recovery_rate(d_hc):
    return tau_max / (1 + d_hc)

# Infection rate based on distance to breeding sites
def infection_rate(d_bs):
    return lambda_base / (1 + d_bs)

# Transmission rate from vectors to humans (adjusted for bed net usage)
def lambda_h(I_v, N_h, t):
    bed_net = bed_net_usage(t)
    return (1 - 0.5*bed_net) * rho_vh * biting_rate( t) * I_v / N_h

# Transmission rate from humans to vectors (consider asymptomatic transmission)
def lambda_v(I_h, A_h, N_h, t):
    bed_net = bed_net_usage(t)
    return (1-0.5*bed_net) * rho_hv * biting_rate( t) * (I_h + 0.5 * A_h) / N_h

# ODE system
def odes(t, y, B, irrigation_factor):
#def odes(t, y):
    S_h, V_h, E_h, A_h, I_h, R_h, S_v, E_v, I_v = y
    N_h = S_h + V_h + E_h + A_h + I_h + R_h  # Total human population

    # Adjusted infection and recovery rates based on proximity
    lam_h = infection_rate(d_bs)*lambda_h(I_v, N_h, t)
    lam_v = infection_rate(d_bs)*lambda_v(I_h, A_h, N_h, t)
    tau_h = recovery_rate(d_hc)  # Recovery rate influenced by proximity to health centers
    p=P(t,A_h,I_h)
    nu_h = nu_h_t(t) 

    # ODE system equations
    dS_h = Lambda_h - lam_h * S_h + alpha_h * R_h + omega_h * V_h - (mu_h + nu_h) * S_h
    dV_h = nu_h * S_h - (1 - 0.75) * lam_h * V_h - (mu_h + omega_h) * V_h
    dE_h = lam_h * S_h + (1 - 0.75) * lam_h * V_h - (gamma_h + mu_h) * E_h
    dA_h = (1 - p) * gamma_h * E_h - (sigma_h + kappa_h + mu_h) * A_h
    dI_h = p * gamma_h * E_h + kappa_h * A_h - (tau_h + delta_h + mu_h) * I_h
    dR_h = tau_h * I_h + sigma_h * A_h - (alpha_h + mu_h) * R_h
    dS_v = r * B * (1 - B / K_b) * irrigation_factor * seasonality(t) - lam_v * S_v - mu_v * S_v
    dE_v = lam_v * S_v - (gamma_v + mu_v) * E_v
    dI_v = gamma_v * E_v - mu_v * I_v

    return [dS_h, dV_h, dE_h, dA_h, dI_h, dR_h, dS_v, dE_v, dI_v]

# Initial conditions
S_h0 = 1880  # Susceptible humans
V_h0 = 0     # Vaccinated humans
E_h0 = 5     # Exposed humans
A_h0 = 2     # Asymptomatic humans
I_h0 = 1     # Infectious humans
R_h0 = 100   # Recovered humans
S_v0 = 4000  # Susceptible vectors
E_v0 = 50    # Exposed vectors
I_v0 = 20    # Infectious vectors

# Time points (in months)
t_span = [0, 1000]   # Simulate for 200 months (approx. 16.67 years)
t_eval = np.linspace(t_span[0], t_span[1], 1000)  # Points to evaluate the solution

###########################################################################################
 #Solve the ODE system
sol = solve_ivp(odes, t_span, [S_h0, V_h0, E_h0, A_h0, I_h0, R_h0, S_v0, E_v0, I_v0], t_eval=t_eval,args=(B, irrigation_factor))
#sol = solve_ivp(odes, t_span, [S_h0, V_h0, E_h0, A_h0, I_h0, R_h0, S_v0, E_v0, I_v0])
# Extract the solution for I_h
I_h_solution = sol.y[4] / (sol.y[0] + sol.y[1] + sol.y[2] + sol.y[3] + sol.y[4] + sol.y[5])
print(len(I_h_solution))
# Plot the infectious human population over time
plt.figure(figsize=(10, 6), dpi=150)
plt.plot(t_eval, I_h_solution, label="Infectious humans (I_h)", color='r')#, marker='o', markersize=4)
plt.xlabel("Time (months)")
plt.ylabel("Prevalence")
#plt.title("Infectious Humans Over Time with Seasonal Variation (in Months)")
#plt.legend()
plt.grid()
plt.show()
##################################################################################################################################################################################
#Define the range of breeding sites and irrigation factors to simulate  
d_hc_range = np.linspace(0, 2, 9)  
d_bs_range = np.linspace(0, 2, 9)  
  
# Initialize a 2D array to store the prevalence values  
prevalence = np.zeros((len(d_hc_range), len(d_bs_range)))  
  
# Simulate the ODE system for each combination of breeding sites and irrigation factors  
for i, d_hc in enumerate(d_hc_range):  
   for j, d_bs in enumerate(d_bs_range):  
      sol = solve_ivp(odes, t_span, [S_h0, V_h0, E_h0, A_h0, I_h0, R_h0, S_v0, E_v0, I_v0], t_eval=t_eval, args=(d_hc, d_bs))  
      I_h_solution = sol.y[4] / (sol.y[0] + sol.y[1] + sol.y[2] + sol.y[3] + sol.y[4] + sol.y[5])  
      prevalence[i, j] = np.mean(I_h_solution)  
  
# Plot the heat map of prevalence  
plt.figure(figsize=(10, 6), dpi=150)
plt.imshow(prevalence, cmap='hot', interpolation='nearest', origin='lower')  
plt.xlabel('distance to breeding site')  
plt.ylabel('distance to health center')
plt.xticks(np.arange(len(d_bs_range)), d_bs_range)  
plt.yticks(np.arange(len(d_hc_range)), d_hc_range)
plt.xticks(rotation=90)
plt.title('Prevalence')  
plt.colorbar()  
plt.show()
########################################################################################
##Define the range of breeding sites and irrigation factors to simulate  
B_range = np.linspace(10, 35, 11)  
irrigation_factor_range = np.linspace(1, 2, 11)  
rirrigation_factor = np.round(irrigation_factor_range, 1)  
# Initialize a 2D array to store the prevalence values  
prevalence = np.zeros((len(B_range), len(irrigation_factor_range)))  
  
# Simulate the ODE system for each combination of breeding sites and irrigation factors  
for i, B in enumerate(B_range):  
   for j, irrigation_factor in enumerate(irrigation_factor_range):  
      sol = solve_ivp(odes, t_span, [S_h0, V_h0, E_h0, A_h0, I_h0, R_h0, S_v0, E_v0, I_v0], t_eval=t_eval, args=(B, irrigation_factor))  
      I_h_solution = sol.y[4] / (sol.y[0] + sol.y[1] + sol.y[2] + sol.y[3] + sol.y[4] + sol.y[5])  
      prevalence[i, j] = np.mean(I_h_solution)  
  
# Plot the heat map of prevalence  
plt.figure(figsize=(10, 6), dpi=150)
plt.imshow(prevalence, cmap='hot', interpolation='nearest', origin='lower')  
plt.xlabel('irrigation factor')  
plt.ylabel('number of breeding sites')
plt.xticks(np.arange(len(rirrigation_factor)), rirrigation_factor)  
plt.yticks(np.arange(len(B_range)), B_range)
plt.xticks(rotation=90)
plt.title('Prevalence')  
plt.colorbar()  
plt.show()



