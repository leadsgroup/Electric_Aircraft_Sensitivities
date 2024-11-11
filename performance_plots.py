import matplotlib.pyplot as plt  
import numpy as np 
    
    
def main():
         
    # -----------------------------------
    P_aircraft   = 0  
    eta_em = 0
    eta_p = 0
    LD = 0
    W0 = 0
    M_struct = 0
    M_pass = 0
    M_crew = 0
    M_payload = 0

    # -----------------------------------    
    n_p = 210
    n_s = 120
    M_cell = 0.048
    E_cell = 0
    e_cell = 0
    Q_cell = 0
    V_cell = 4.2
    alpha = 0

    # -----------------------------------    
    P_motor = 0
    T_torque_density = 0
    omega = 0
    Pd_motor_cooling = 0
    eta_motor = 0

    # -----------------------------------    
    P_inverter = 0
    Pd_inverter = 0
    eta_inverter = 0

    # -----------------------------------    
    eta_propeller = 0.86
    eta_battery = 0.9
    Pd_inverter_cooling = 0
    V = 0
    E0 = 2.51e3
    r_cond = 0
    rho = 0
    rho_theta_insul = 0
    L = 0
    rho_cond = 0
    rho_insul = 0
     
    # -----------------------------------    
    theta_a = 253 # in Kelvin
    I = 0
    T_4 = 0
    D = 0
    g = 0
    M_0 = 0
             
    Range,E_pack =  compute_performance(n_p, n_s, M_cell, E_cell, e_cell, alpha, P_motor, T_torque_density, omega, Pd_motor_cooling, eta_motor, 
             P_inverter, Pd_inverter, eta_inverter, eta_propeller, eta_battery, P_aircraft, Pd_inverter_cooling, V, E0, r_cond, rho, rho_theta_insul, 
             L, rho_cond, rho_insul, eta_em, eta_p, LD, W0, M_struct, M_pass, M_crew, M_payload, theta_a, I, T_4, D, g, M_0, Q_cell, V_cell)            
 

    return

def compute_performance(n_p, n_s, M_cell, E_cell, e_cell, alpha, P_motor, T_torque_density, omega, Pd_motor_cooling, eta_motor, 
         P_inverter, Pd_inverter, eta_inverter, eta_propeller, eta_battery, P_aircraft, Pd_inverter_cooling, V, E0, r_cond, rho, rho_theta_insul, 
         L, rho_cond, rho_insul, eta_em, eta_p, LD, W0, M_struct, M_pass, M_crew, M_payload, theta_a, I, T_4, D, g, M_0, Q_cell, V_cell):

     
    # SECTION I: Electrochemical Energy Storage Systems (Batteries)
    M_energy_storage =  energy_storage()

    # SECTION II: Electrical Power Conversion Machines
    M_electric_power_conversion =  power_conversion()
    
    # SECTION III: Cabling Mass
    M_cable =  cable_mass()

    # SECTION IV: Overall Performance

    # Equation (25): Total Drivetrain Mass (M_drivetrain) is the sum of storage, power conversion, and cabling masses
    M_drivetrain = M_energy_storage + M_electric_power_conversion + M_cable # Equation (25)

    # Equation (24): Total Aircraft Mass (M_0) including masses of structure, passengers, crew, payload, and drivetrain
    M_0 = M_struct + M_pass + M_crew + M_payload + M_drivetrain # Equation (24)

    # Equation (23): Aircraft Range Equation
    Range =  aircraft_flight_range()
    
    # Equation (27a): Total Battery Pack Charge
    Q_pack = n_p * Q_cell  # Equation (27a)

    # Equation (27b): Total Battery Pack Voltage (V_pack)
    V_pack = n_s* V_cell  # Equation (27b)

    # Equation (26): Total Battery Pack Energy (E_pack), derived from charge and voltage
    E_pack = Q_pack * V_pack # Equation (26)

    return Range,  E_pack

def energy_storage(n_p, n_s, M_cell, E_cell, e_cell, alpha, P_motor, T_torque_density, omega, Pd_motor_cooling, eta_motor, 
         P_inverter, Pd_inverter, eta_inverter, eta_propeller, eta_battery, P_aircraft, Pd_inverter_cooling, V, E0, r_cond, rho, rho_theta_insul, 
         L, rho_cond, rho_insul, eta_em, eta_p, LD, W0, M_struct, M_pass, M_crew, M_payload, theta_a, I, T_4, D, g, M_0, Q_cell, V_cell):
    
    # Equation (2): Total dry battery mass (higher fidelity for Eq. (1))
    M_pack = (E_cell / e_cell) * n_p * n_s  # Equation (2)

    # Equation (3): Battery Control & Thermal Management System Mass
    M_BMS = alpha * M_pack  # Equation (3)

    # Equation (4): Total Energy Storage Mass
    M_energy_storage = (1 + alpha) * (E_cell / e_cell) * n_p * n_s
    
    return M_energy_storage

def power_conversion(n_p, n_s, M_cell, E_cell, e_cell, alpha, P_motor, T_torque_density, omega, Pd_motor_cooling, eta_motor, 
         P_inverter, Pd_inverter, eta_inverter, eta_propeller, eta_battery, P_aircraft, Pd_inverter_cooling, V, E0, r_cond, rho, rho_theta_insul, 
         L, rho_cond, rho_insul, eta_em, eta_p, LD, W0, M_struct, M_pass, M_crew, M_payload, theta_a, I, T_4, D, g, M_0, Q_cell, V_cell): 

    # Equation (6): Motor Mass, with motor power calculated based on aircraft power, motor efficiency, and propeller efficiency
    P_motor = P_aircraft / (eta_motor * eta_propeller)
    M_motor = P_motor / (T_torque_density * omega)  # Equation (6)

    # Equation (7) or (9): Motor Cooling System Mass (using higher fidelity if needed)
    P_inverter = P_aircraft / (eta_inverter * eta_motor * eta_propeller)  # Derived power of the inverter based on efficiencies
    M_motor_cooling = P_inverter * (1 - eta_motor) / Pd_motor_cooling  # Equation (7) or Equation (9) as per fidelity

    # Equation (10): Inverter Mass
    M_inverter = P_inverter / Pd_inverter  # Equation (10)
    # Alternative: M_inverter = (k_1 / f_s) + k_2 * f_s * V**2 for different fidelity level (Eq. (13))

    # Equation (11) or (12): Inverter Cooling System Mass (using higher fidelity if needed)
    P_battery = P_aircraft / (eta_battery * eta_inverter * eta_motor * eta_propeller)  # Battery power derived based on efficiencies
    M_inverter_cooling = P_battery * (1 - eta_inverter) / Pd_inverter_cooling  # Equation (11) or Equation (12)

    # Equation (15): Total Motor Drive Mass
    M_total_motor_drive = M_motor + M_motor_cooling + M_inverter + M_inverter_cooling  # Equation (15)

    # Equation (16): Electric Power Conversion Mass (using detailed formula for high fidelity)
    M_electric_power_conversion = (P_aircraft / (eta_motor * eta_propeller)) * ((1 / T_torque_density) + ((1-eta_motor) / (eta_inverter * Pd_motor_cooling)) + (1 / (eta_inverter * Pd_inverter)) + ((1 - eta_inverter) / (eta_battery * eta_inverter * Pd_inverter_cooling)))  # Equation (16)

    return M_electric_power_conversion  

def cable_mass(n_p, n_s, M_cell, E_cell, e_cell, alpha, P_motor, T_torque_density, omega, Pd_motor_cooling, eta_motor, 
         P_inverter, Pd_inverter, eta_inverter, eta_propeller, eta_battery, P_aircraft, Pd_inverter_cooling, V, E0, r_cond, rho, rho_theta_insul, 
         L, rho_cond, rho_insul, eta_em, eta_p, LD, W0, M_struct, M_pass, M_crew, M_payload, theta_a, I, T_4, D, g, M_0, Q_cell, V_cell):

    # Equation (18): Cable Insulation Radius based on voltage and electric field constraints
    # E0 is the electric field
    r_insul = r_cond * np.exp(V / (E0 * r_cond))  # Equation (18)

    # Equation (20): Conductor Resistance (thermal constraint based on material properties)
    R_prime = rho / (np.pi * r_cond ** 2)  # Equation (20)

    # Equation (21): Thermal Resistance of the insulation
    T_1 = rho_theta_insul / (2 * np.pi) * np.log(r_insul / r_cond)  # Equation (21)

    # Equation (22): Total Cable Mass calculation based on conductor and insulation volume and density
    M_cable = np.pi * L * (r_cond ** 2 * rho_cond + (r_insul ** 2 - r_cond ** 2) * rho_insul)  # Equation (22)

    # Equation (19): Maximum Temperature (conductor temperature based on current, resistance, and thermal resistances)
    theta_max = theta_a + I**2 * R_prime * (T_1 + T_4)  # Equation (19)
    
    return M_cable

def aircraft_flight_range(n_p, n_s, M_cell, E_cell, e_cell, alpha, P_motor, T_torque_density, omega, Pd_motor_cooling, eta_motor, 
         P_inverter, Pd_inverter, eta_inverter, eta_propeller, eta_battery, P_aircraft, Pd_inverter_cooling, V, E0, r_cond, rho, rho_theta_insul, 
         L, rho_cond, rho_insul, eta_em, eta_p, LD, W0, M_struct, M_pass, M_crew, M_payload, theta_a, I, T_4, D, g, M_0, Q_cell, V_cell):

    R = eta_em * eta_p * (L/D) * E_pack / (g * M_0)  # Equation (23)
     
    return R 

if __name__ == '__main__':
    main() 
    plt.show()    