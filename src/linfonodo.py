import numpy as np

def diferential(y, parameters):
    # DT Dendritic cells in the tissue
    # FL Average number of antibodies on tissue 
    V_LV = parameters[0]     # Volume of the domain areas in contact with lymph vessels
    V_LN = parameters[1]     # Volume of the lymph node
    V_BV = parameters[2]
    
    gamma_D = parameters[3] # Dendritic cells migration rate (Tissue -> Lymph node)
    alpha_T_c = parameters[4]  
    estable_T_c = parameters[5]
    gamma_T = parameters[6]    
    theta_BV = parameters[7]
    Tt_c = parameters[8] # T citotoxic from tissue

    b_T = parameters[9]    
    rho_T = parameters[10]    
    b_rho = parameters[11]    
    alpha_T_h = parameters[12]    
    estable_T_h = parameters[13]
    rho_B = parameters[14]    #* Replication rate of B cells
    alpha_B = parameters[15]  #* Replication rate of B cells 
    estable_B = parameters[16]# Estable B value
    
    rho_F = parameters[17]  # Antibodies production rate (by B cells) 
    gamma_F = parameters[18]# Antibodies migration rate (to tissue)
    DT = parameters[19]
    FT = parameters[20]
    b_rho_b = parameters[21]

    # When working with only one diferential function there must be created a null extra position
    dy = np.zeros(5)

    # Dendritic cells
    dy[0] = gamma_D * (DT - y[0]) * (V_LV / V_LN)
    # Cytotoxic T cells
    dy[1] = alpha_T_c * (estable_T_c - y[1]) - (gamma_T * theta_BV * (y[1] - Tt_c)) * (V_BV / V_LN)
    # Helper T cells
    dy[2] = b_T*(rho_T * y[2] * y[0] - y[2]*y[0]) - (b_rho * y[2] * y[0] * y[3]) + alpha_T_h * (estable_T_h - y[2])
    # B cells
    dy[3] = (b_rho_b * ((rho_B * y[2] * y[0]) - (y[2] * y[0] * y[3]))) + alpha_B * (estable_B - y[3])
    # Antibodies
    dy[4] = rho_F * y[3] - ((gamma_F * theta_BV * (y[4] - FT)) * (V_BV / V_LN))

    return dy