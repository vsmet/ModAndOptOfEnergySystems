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
#param temp_supply_high{b in BUILDINGS,t in TIME} >= 0; #deg C
#param temp_return_high{b in BUILDINGS,t in TIME} >= 0; #deg C

/*******************************************************/
# Demand parameters
/*******************************************************/
param spec_annual_heat_demand{b in BUILDINGS} >= 0, default 0;    #kJ/m2(yr)
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
subject to Heat_Demand_Constr{b in BUILDINGS, t in TIME}:
    Heat_Demand[b,t] = 
    if (external_temp[t] < temp_threshold[b]) then
      (k1[b] * (external_temp[t]) + k2[b])*1000            #kW
    else 0;
  

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

param  f_act   := Year_ind/Ref_ind;

var Capacity {c in COMPONENTS} >= 0;
var PC{c in COMPONENTS} >= 0;
var BM_C{c in COMPONENTS} >= 0;
var GR_C{c in COMPONENTS} >= 0;
var an_CAPEX {c in COMPONENTS} >=0;
var an_CAPEX_Tot  >=0;


############################################################################################
# CONSTRAINTS
############################################################################################

# BOILER MODEL
param Efficiency_Boiler := 0.98;


var NG_Demand_Boiler{t in TIME} >= 0;
var Heat_Supple_Boiler{t in TIME} >= 0;

#Energy model
subject to Boiler_Energy_Balance_Constr{t in TIME}:
  Heat_Supple_Boiler[t] = Efficiency_Boiler*NG_Demand_Boiler[t];  #kW  
subject to Boiler_Size_Constr{t in TIME}:
  Heat_Supple_Boiler[t] <= Capacity["BOILER"];             #kW

#*********** HP MODEL *****************

param COP{b in BUILDINGS, t in TIME} := (7.2-(7.2-4.7)/(20*(-30+temp_supply[b,t])));

var EL_Demand_HP{b in BUILDINGS,t in TIME} >=0;
var Heat_Supple_HP{b in BUILDINGS,t in TIME} >= 0;
var Capacity_HP{b in BUILDINGS} >= 0;

#Energy model
subject to HP_Energy_Balance_Constr{b in BUILDINGS,t in TIME}:
  Heat_Supple_HP[b,t] = COP[b,t]*EL_Demand_HP[b,t];  #kW  
subject to HP_Size_Constr{b in BUILDINGS,t in TIME}:
  Heat_Supple_HP[b,t] <= Capacity_HP[b];             #kW
subject to HP_Size_1:
  Capacity["HEATPUMPHIGH"] = Capacity_HP["EPFLhigh"];
subject to HP_Size_2:
  Capacity["HEATPUMPLOW"] = Capacity_HP["EPFLlow"];  


#SOLAR PANEL MODEL#########################################
param Efficiency_SolarPanels := 0.11327;  #Voir feuille excel DATA, Sylvain
param solarfarm_area := 15500;            #m^2, donnée du projet

var El_Available_Solar{t in TIME} >=0;

subject to El_available_Constr{t in TIME}:
  El_Available_Solar[t] = solarfarm_area*Efficiency_SolarPanels*solar_radiation[t]; #kW
  

#MASS BALANCE NATURAL GAS
var NG_Demand_grid{t in TIME} >= 0;

subject to Natural_gas_Demand_Constr{t in TIME}:
  NG_Demand_grid[t] = NG_Demand_Boiler[t];    #kW
  

# HEAT BALANCE 
subject to MidHigh_Circuit_Constr{b in BUILDINGS,t in TIME}:
  Heat_Supple_HP[b,t] <= Heat_Demand[b,t];
subject to Heat_balance_Constr{t in TIME}:
  Heat_Supple_Boiler[t]= sum{b in BUILDINGS}( Heat_Demand[b,t] - Heat_Supple_HP[b,t]);   #kW
subject to EightyeightPerc_Constr:
  sum{b in BUILDINGS, t in TIME}(Heat_Supple_HP[b,t]) >= 0.88*sum{b in BUILDINGS,t in TIME}(Heat_Demand[b,t]); #SYSTEM REQUIREMENTS
subject to Peak:
  25000 - Capacity["HEATPUMPLOW"] - Capacity["HEATPUMPHIGH"] <= Capacity["BOILER"]; 



#ELECTRICITY BALANCE
var El_Buy{t in TIME} >=0;
var El_Sell{t in TIME}>=0;
subject to Electricity_balance_Constr{t in TIME}:
  El_Available_Solar[t] + El_Buy[t] - El_Sell[t] - sum{b in BUILDINGS} EL_Demand_HP[b,t]= sum{b in BUILDINGS} Elec_Demand[b,t]; #kW

#INVESTMENT
subject to PC_Con{c in COMPONENTS}:
  PC[c] = f_act*((PC_max[c] - PC_min[c])*((Capacity[c] - C_min[c])/(C_max[c] - C_min[c])) + PC_min[c]); #USD
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
display opex;
display Capacity;


end;
