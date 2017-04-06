############################################################################################
# Course: Modelling and optimisation of energy systems
# course spring semester 2017
# EPFL Campus optimization

# IPESE, EPFL

############################################################################################
# SETS
############################################################################################
set HP;
set TIME;
set COMPONENTS;
set HEAT_UTIL;
set ELEC_UTIL;
set FUEL_USERS;
set NG_USERS;



############################################################################################
# PARAMETER
############################################################################################

/*******************************************************/
# General parameters
/*******************************************************/
param TIMEsteps{t in TIME};         #hr

/*******************************************************/
# Meteo parameters
/*******************************************************/
param external_temp{t in TIME};       #deg C
param solar_radiation{t in TIME};     # kw/m2

/*******************************************************/
# Building parameters
/*******************************************************/
param floor_area;        #m2
param temp_threshold;       #deg C
param temp_supply{b in HP,t in TIME} >= 0; #deg C
param temp_return{b in HP,t in TIME} >= 0; #deg C

/*******************************************************/
# Demand parameters
/*******************************************************/
param spec_annual_elec_demand;    #kWh/m2(yr)






############################################################################################
# VARIABLES (and defining equations)
############################################################################################
# Building model using 
# - area specific energy demand data to determine building demand
# - Energy Signature (ES) to determine power demand


# COOLING LOOP VARIABLES
param vol_cooling_water{t in TIME}; # m3/mois
param pumping_cost := 0.304 ; #kWh/m3,  calculé depuis http://exploitation-energies.epfl.ch/bilans_energetiques/details/cct, Sylvain
param cp_water:=4.18; #kJ/kg K
# Parameter heating signature
param k1{b in HP};
param k2{b in HP};
 
# TIME-DEPENDENT HEAT DEMAND
var Heat_Demand{b in HP, t in TIME} >= 0;
subject to Heat_Demand_Constr{b in HP, t in 1..12}:
    Heat_Demand[b,t] = 
    if (external_temp[t] < temp_threshold) then
      (k1[b] * (external_temp[t]) + k2[b])*1000            #kW
    else 0;

# INVESTMENT VARIABLES
param ex_USD_CHF;
param interest_rate;
param lifetime;
param Year_ind;
param Ref_ind;
param F_P{c in COMPONENTS};
param F_T{c in COMPONENTS};
param F_BM {c in COMPONENTS};
param alpha_1;
param alpha_2;
param C_min{c in COMPONENTS};
param C_max{c in COMPONENTS};
param PC_min {c in COMPONENTS};
param PC_max {c in COMPONENTS};  
param f_act := Year_ind/Ref_ind;

var Capacity {c in COMPONENTS} >= 0;
var PC{c in COMPONENTS} >= 0;
var BM_C{c in COMPONENTS} >= 0;
var GR_C{c in COMPONENTS} >= 0;
var an_CAPEX {c in COMPONENTS} >=0;
var an_CAPEX_Tot  >=0;

# COMPONENT VARIABLES
param Component_temp {c in COMPONENTS, t in TIME};
var Component_Use {c in COMPONENTS} binary;
var ComponentSize_t {c in COMPONENTS, t in TIME} >= 0;

subject to Size_Constr{c in COMPONENTS, t in TIME}:
  ComponentSize_t[c,t] <= Capacity[c];             #kW

subject to Component_cmin_cstr {c in COMPONENTS}:
  Component_Use[c]*C_min[c] <= Capacity[c];

subject to Component_cmax_cstr {c in COMPONENTS}:
  Component_Use[c]*C_max[c] >= Capacity[c];

############################################################################################
# CONSTRAINTS
############################################################################################

# HEAT BALANCE ########################################
#Energy model
var Heating_LT {c in HEAT_UTIL, t in TIME} >= 0;
var Heating_HT {c in HEAT_UTIL, t in TIME} >= 0;

subject to Heat_Demand_2050{b in HP}:
  Heat_Demand[b,13] = C_max[b];
subject to Energy_Balance_LT_cstr2 {t in TIME, bu in HP: temp_supply[bu,t]<=50}: 
  sum{b in HP} Heat_Demand[b,t] = sum {c in HEAT_UTIL} ComponentSize_t [c,t];
# Energy balance for LT HP
subject to Energy_Balance_LT_cstr {b in HP,t in TIME: temp_supply[b,t]>50}:
  Heat_Demand['HPLOW',t] = sum {c in HEAT_UTIL: c!="HPHIGH"} Heating_LT [c,t];
# Energy balance for HT HP
subject to Energy_Balance_HT_cstr {b in HP,t in TIME: temp_supply[b,t]>50}:
  Heat_Demand['HPHIGH',t] = sum {c in HEAT_UTIL : c!="HPLOW"} Heating_HT [c,t];


# Overall energy balance
subject to Energy_Balance_overall_cstr {c in HEAT_UTIL,t in TIME}:
  ComponentSize_t[c,t] = Heating_LT[c,t]+Heating_HT[c,t];

#88% constraint
subject to EightyeightPerc_Constr:
  sum{h in HP, t in TIME}(ComponentSize_t[h,t]*TIMEsteps[t]) >= 0.88*sum{b in HP,t in TIME}(Heat_Demand[b,t]*TIMEsteps[t]); #SYSTEM REQUIREMENTS

# COP and FUEL_using efficiencies
param lake_temp := 7;
param carnot_eff := 0.5;
param COP_th{h in HP, t in TIME} := (Component_temp[h,t]+273)/(Component_temp[h,t]-lake_temp);
param COP{h in HP, t in TIME} := carnot_eff * COP_th [h,t];

param FUEL_el_eff{u in FUEL_USERS}; 
param FUEL_th_eff{u in FUEL_USERS}; 

 
# FUEL BALANCE  #################################################
var FUEL_Demand{u in FUEL_USERS, t in TIME}>=0;

subject to FUEL_heat_balance_constr{u in FUEL_USERS, t in TIME}: 
 ComponentSize_t[u,t] = FUEL_th_eff[u]*FUEL_Demand[u,t]; #dim

# ELEC BALANCE #################################################
var El_prod{u in COMPONENTS, t in TIME};
var El_Buy{t in TIME} >=0;
var El_Sell{t in TIME}>=0;
var Elec_Demand{t in TIME}>=0;

param Efficiency_SolarPanels := 0.11327;  #Voir feuille excel DATA, Sylvain
param solarfarm_area := 15500;            #m^2, donnée du projet

subject to Non_elec_prod_constr{t in TIME}:
  El_prod["BOILER",t] = 0;
subject to FUEL_elec_balance_constr{u in FUEL_USERS, t in TIME}: 
  El_prod[u,t] = FUEL_el_eff[u]*FUEL_Demand[u,t]; #dim
subject to El_available_Constr{t in TIME}: #AJOUTER COUTS DES NOUVEAUX PANNEAUX!
  El_prod["SOLAR",t] = ((solarfarm_area+Capacity["SOLAR"])*Efficiency_SolarPanels)*solar_radiation[t]; #kW
subject to HP_Energy_Balance_cstr{h in HP,t in TIME}:
  ComponentSize_t[h,t] = COP[h,t]*(-El_prod[h,t]);  #kW
subject to Elec_demand_system{t in TIME}:
Elec_Demand[t] = (spec_annual_elec_demand*floor_area)/12 + (vol_cooling_water[t]+(Heat_Demand['HPLOW',t]+Heat_Demand['HPHIGH',t])/(4*cp_water))*pumping_cost; #dim

subject to Electricity_balance_Constr{t in TIME}:
  sum{u in COMPONENTS} El_prod[u,t]+ El_Buy[t] - El_Sell[t]= Elec_Demand[t]; #kW



# INVESTMENT ###################################################
subject to PC_Con{c in COMPONENTS}:
  PC[c] = f_act*((PC_max[c] - PC_min[c])*(Capacity[c]/(C_max[c] - C_min[c])) + Component_Use[c]*PC_min[c]); #USD
subject to BM_C_Con{c in COMPONENTS}:
  BM_C[c] = F_P[c]*F_T[c]*F_BM[c]*PC[c];
subject to GR_C_Con{c in COMPONENTS}:
  GR_C[c] = BM_C[c]*(alpha_1*alpha_2 + 1);    #Cout totaux
subject to an_CAPEX_Con{c in COMPONENTS}:
  an_CAPEX[c] = GR_C[c]*((interest_rate*(1+interest_rate)^lifetime)/((1+interest_rate)^lifetime - 1));  #Cout annualisés
subject to an_CAPEXTot_Con:
  an_CAPEX_Tot = sum{c in COMPONENTS} an_CAPEX[c];



  
############################################################################################
# OBJECTIVE FUNCTION
############################################################################################

/*******************************************************/
# Economic parameters
/*******************************************************/
param c_el_in;
param c_el_out;
param c_ng_in;
param c_ds_in;


minimize opex:
sum{u in NG_USERS,t in TIME} ((c_ng_in*FUEL_Demand[u,t] + c_ds_in*FUEL_Demand["ICENGINE",t] + c_el_in*El_Buy[t] - c_el_out*El_Sell[t])*TIMEsteps[t]) + an_CAPEX_Tot;

solve;

# To do!
display Heat_Demand;
display opex;
display Capacity;
display Component_Use;
display Heat_Demand;

display ComponentSize_t;

display Capacity;
display an_CAPEX;


display El_Buy;

#display COP;

end;
