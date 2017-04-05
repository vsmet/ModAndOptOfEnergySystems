############################################################################################
# Course: Modelling and optimisation of energy systems
# course spring semester 2017
# EPFL Campus optimization

# IPESE, EPFL

############################################################################################
# SETS
############################################################################################
set BUILDINGS;
set TIME;
set COMPONENTS;


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
param floor_area{b in BUILDINGS} >= 0;        #m2
param temp_threshold{b in BUILDINGS};       #deg C
param temp_supply{b in BUILDINGS,t in TIME} >= 0; #deg C
param temp_return{b in BUILDINGS,t in TIME} >= 0; #deg C

/*******************************************************/
# Demand parameters
/*******************************************************/
#param spec_annual_heat_demand{b in BUILDINGS} >= 0, default 0;    #kJ/m2(yr)
param spec_annual_elec_demand{b in BUILDINGS} >= 0, default 0;    #kWh/m2(yr)

############################################################################################
# VARIABLES (and defining equations)
############################################################################################
# Building model using 
# - area specific energy demand data to determine building demand
# - Energy Signature (ES) to determine power demand

/*******************************************************/
# Energy variables
/*******************************************************/ 
# ELEC
var Annual_Elec_Demand{b in BUILDINGS} >= 0;
subject to Annual_Elec_Demand_Constr{b in BUILDINGS}:
  Annual_Elec_Demand[b] = floor_area[b] * spec_annual_elec_demand[b];   #kWh(/yr)
  
# Parameter heating signature
param k1{b in BUILDINGS};
param k2{b in BUILDINGS};
 

# TIME-DEPENDENT HEAT DEMAND
var Heat_Demand{b in BUILDINGS, t in TIME} >= 0;
subject to Heat_Demand_Constr{b in BUILDINGS, t in 1..12}:
    Heat_Demand[b,t] = 
    if (external_temp[t] < temp_threshold[b]) then
      (k1[b] * (external_temp[t]) + k2[b])*1000            #kW
    else 0;

param max_demand{b in BUILDINGS};
subject to Heat_Demand_2050{b in BUILDINGS}:
  Heat_Demand[b,13] = max_demand[b];



# TIME-DEPENDENT ELEC DEMAND
var Elec_Demand{b in BUILDINGS, t in TIME} >= 0;
subject to Elec_Demand_Constr{b in BUILDINGS, t in TIME}:
  Elec_Demand[b,t] = Annual_Elec_Demand[b] / 12;    #kW

/*******************************************************/
# Investment variables
/*******************************************************/ 
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

# ENERGY BALANCE ########################################
#Energy model
var Heating_LT {c in COMPONENTS, t in TIME} >= 0;
var Heating_HT {c in COMPONENTS, t in TIME} >= 0;

# Energy balance for LT buildings
subject to Energy_Balance_LT_cstr {b in BUILDINGS,t in TIME}:
  Heat_Demand['EPFLlow',t] = sum {c in COMPONENTS} (if Component_temp[c,t]>= (temp_supply[b,t]+3) then (Heating_LT [c,t]) else (0));

# Energy balance for HT buildings
subject to Energy_Balance_HT_cstr {b in BUILDINGS,t in TIME}:
  Heat_Demand['EPFLhigh',t] = sum {c in COMPONENTS} (if Component_temp[c,t]>= (temp_supply[b,t] + 5) then (Heating_HT [c,t]) else (0));

# Overall energy balance
subject to Energy_Balance_overall_cstr {c in COMPONENTS,t in TIME}:
  ComponentSize_t[c,t] = Heating_LT[c,t]+Heating_HT[c,t];


# HP MODEL ################################################
set HP;
param lake_temp := 7;
param carnot_eff := 0.5;
param COP_th{h in HP, t in TIME} := (Component_temp[h,t]+273)/(Component_temp[h,t]-lake_temp);
param COP{h in HP, t in TIME} := carnot_eff * COP_th [h,t];

# Energy balance for HP
var EL_Demand_HP {h in HP, t in TIME};
subject to HP_Energy_Balance_cstr{h in HP,t in TIME}:
  ComponentSize_t[h,t] = COP[h,t]*EL_Demand_HP[h,t];  #kW 

/*# SOFC MODEL ##############################################
set COGENERATION;
param SOFC_el_eff := 0.2; #dimitri inventé
param SOFC_th_eff := 0.5; #dimitri inventé

# Energy balance for SOFC
var NG_Demand_SOFC{t in TIME}>=0;
var El_prod_SOFC{t in TIME} >=0;

subject to SOFC_heat_balance_constr{t in TIME}: 
 ComponentSize_t['SOFC',t] = SOFC_th_eff*NG_Demand_SOFC[t]; #dim
subject to SOFC_elec_balance_constr{t in TIME}: 
 El_prod_SOFC[t] = SOFC_el_eff*NG_Demand_SOFC[t]; #dim*/

 #§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§

set COGENERATION;
param COG_el_eff{u in COGENERATION}; 
param COG_th_eff{u in COGENERATION}; 

# Energy balance for SOFC
var NG_Demand_COG{u in COGENERATION, t in TIME}>=0;
var El_prod_COG{u in COGENERATION, t in TIME} >=0;

subject to COG_heat_balance_constr{u in COGENERATION, t in TIME}: 
 ComponentSize_t[u,t] = COG_th_eff[u]*NG_Demand_COG[u,t]; #dim
subject to COG_elec_balance_constr{u in COGENERATION, t in TIME}: 
 El_prod_COG[u,t] = COG_el_eff[u]*NG_Demand_COG[u,t]; #dim

# BOILER MODEL  ###########################################
param Efficiency_Boiler := 0.98;
var NG_Demand_Boiler{t in TIME} >= 0;

# Energy balance for Boiler
subject to Boiler_Energy_Balance_Constr{t in TIME}:
  ComponentSize_t['BOILER',t] = Efficiency_Boiler*NG_Demand_Boiler[t];  #kW  



# SOLAR PANEL MODEL#########################################
param Efficiency_SolarPanels := 0.11327;  #Voir feuille excel DATA, Sylvain
param solarfarm_area := 15500;            #m^2, donnée du projet
var solarfarm_area_increase >= 0, <= 3500; #m^2, Valeur max estimated quickly from Maps
var El_Available_Solar{t in TIME} >=0;

subject to El_available_Constr{t in TIME}: #AJOUTER COUTS DES NOUVEAUX PANNEAUX!
  El_Available_Solar[t] = (solarfarm_area*Efficiency_SolarPanels+solarfarm_area_increase*Efficiency_SolarPanels)*solar_radiation[t]; #kW

# MASS BALANCE NATURAL GAS #################################
var NG_Demand_grid{t in TIME} >= 0;

subject to Natural_gas_Demand_Constr{t in TIME}:
  NG_Demand_grid[t] = NG_Demand_Boiler[t]+sum{u in COGENERATION} NG_Demand_COG[u,t];  #kW #dimitri

# COOLING LOOP MODEL ########################################
param vol_cooling_water{t in TIME}; # m3/mois
param pumping_cost := 0.304 ; #kWh/m3,  calculé depuis http://exploitation-energies.epfl.ch/bilans_energetiques/details/cct, Sylvain
#param cooling_water_Tin = 6 ; # °C, 1.7.1. page 5
#param cooling_water_Tout = 13; # °C, 1.7.1. page 5
var El_pump_cooling{t in TIME} >= 0;

subject to El_pump_cooling_water{t in TIME}:
  El_pump_cooling[t]=vol_cooling_water[t]*pumping_cost;
 

# HEAT BALANCE #################################################
subject to EightyeightPerc_Constr:
  sum{h in HP, t in TIME}(ComponentSize_t[h,t]) = 0.88*sum{b in BUILDINGS,t in TIME}(Heat_Demand[b,t]); #SYSTEM REQUIREMENTS
# subject to Peak:
#   25000 - Capacity["HEATPUMPLOW"] - Capacity["HEATPUMPHIGH"] <= Capacity["BOILER"]; 


# ELECTRICITY BALANCE ##########################################
var El_Buy{t in TIME} >=0;
var El_Sell{t in TIME}>=0;
subject to Electricity_balance_Constr{t in TIME}:
  El_Available_Solar[t]+ sum{u in COGENERATION} El_prod_COG[u,t]+ El_Buy[t] - El_Sell[t] - sum{h in HP} EL_Demand_HP[h,t] - El_pump_cooling[t]= sum{b in BUILDINGS} Elec_Demand[b,t]; #kW

# INVESTMENT ###################################################
subject to PC_Con{c in COMPONENTS}:
  PC[c] = f_act*((PC_max[c] - PC_min[c])*(Capacity[c]/(C_max[c] - C_min[c])) + Component_Use[c]*PC_min[c]); #USD
subject to BM_C_Con{c in COMPONENTS}:
  BM_C[c] = F_P[c]*F_T[c]*F_BM[c]*PC[c];
subject to GR_C_Con{c in COMPONENTS}:
  GR_C[c] = BM_C[c]*(alpha_1*alpha_2 + 1);
subject to an_CAPEX_Con{c in COMPONENTS}:
  an_CAPEX[c] = GR_C[c]*((interest_rate*(1+interest_rate)^lifetime)/((1+interest_rate)^lifetime - 1));
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


minimize opex:
sum{t in TIME}((c_ng_in*NG_Demand_grid[t] + c_el_in*El_Buy[t] - c_el_out*El_Sell[t])*TIMEsteps[t]) + an_CAPEX_Tot;



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

display COP;

end;
