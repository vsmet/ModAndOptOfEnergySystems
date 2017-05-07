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

param TIMEsteps{t in TIME};         #hr
param external_temp{t in TIME};       #deg C
param solar_radiation;     # kw/m2
param floor_area;        #m2
param temp_threshold;       #deg C
param spec_annual_elec_demand;    #kWh/m2(yr)
param vol_cooling_water; # m3/mois
param pumping_cost := 0.304 ; #kWh/m3,  calculé depuis http://exploitation-energies.epfl.ch/bilans_energetiques/details/cct, Sylvain
param cp_water:=4.18; #kJ/kg K
# Parameter heating signature
param k1{b in HP};
param k2{b in HP};
param temp_supply {h in HP, t in TIME} := 
  if h = "HPHIGH"
    then (-1.431*external_temp[t] + 50.769) else (-0.9231*external_temp[t] + 40.769); 
param HPTemp {h in HP, t in TIME} :=temp_supply[h,t] + 5;
param Heat_Demand { h in HP, t in TIME} :=
  if (external_temp[t] < temp_threshold) then
      (k1[h] * (external_temp[t]) + k2[h])*1000            #kW
    else 0;
# Parameter Investment
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
param Component_temp {c in COMPONENTS, t in TIME};
# COP and FUEL_using efficiencies
param lake_temp;
param carnot_eff := 0.5;
param COP_th{h in HP, t in TIME} := (HPTemp[h,t]+273)/(HPTemp[h,t]-lake_temp);
param COP{h in HP, t in TIME} := carnot_eff * COP_th [h,t];
param FUEL_el_eff{u in FUEL_USERS}; 
param FUEL_th_eff{u in FUEL_USERS}; 
#Solar Farm
param Efficiency_SolarPanels;
param solarfarm_area;
#Cost
param c_el_in;
param c_el_out;
param c_ng_in;
param c_ds_in;

############################################################################################
# VARIABLES (and defining equations)
############################################################################################

#Variables Investment
var Capacity {c in COMPONENTS} >= 0;
var PC{c in COMPONENTS} >= 0;
var BM_C{c in COMPONENTS} >= 0;
var GR_C{c in COMPONENTS} >= 0;
var an_CAPEX {c in COMPONENTS} >=0;
var an_CAPEX_Tot  >=0;
var Component_Use {c in COMPONENTS} binary;
var ComponentSize_t {c in COMPONENTS, t in TIME} >= 0;
var Heating_LT {c in HEAT_UTIL, t in TIME} >= 0;
var Heating_HT {c in HEAT_UTIL, t in TIME} >= 0;
var FUEL_Demand{u in FUEL_USERS, t in TIME}>=0;
#Variables Electricity
var El_prod{u in COMPONENTS, t in TIME};
var El_Buy{t in TIME} >=0;
var El_Sell{t in TIME}>=0;
var Elec_Demand{t in TIME} >=0;

############################################################################################
# CONSTRAINTS
############################################################################################

#Component Sizing
subject to Size_Constr{c in COMPONENTS, t in TIME}:
  ComponentSize_t[c,t] <= Capacity[c];             #kW
subject to Component_cmin_cstr {c in COMPONENTS}:
  Component_Use[c]*C_min[c] <= Capacity[c];
subject to Component_cmax_cstr {c in 
  COMPONENTS}:
    Component_Use[c]*C_max[c] >= Capacity[c];

#Constraints For temperatures of HeatpumpCycles
subject to Energy_Balance_LT_cstr1 {t in TIME: temp_supply['HPLOW',t]<=50}: 
  sum{b in HP} Heat_Demand[b,t] = sum {c in HEAT_UTIL} ComponentSize_t [c,t];

subject to Energy_Balance_LT_cstr2 { t in TIME: temp_supply['HPLOW',t]>50}:
  Heat_Demand['HPLOW',t] = sum {c in HEAT_UTIL: c!= 'HPLOW'} Heating_LT [c,t];

subject to Energy_Balance_LT_cstr4 { t in TIME: temp_supply['HPHIGH',t] < 65 && temp_supply['HPHIGH',t]>50}:
  Heat_Demand['HPHIGH',t] = sum {c in HEAT_UTIL: c!= 'HPLOW'} Heating_HT [c,t];

#/!\ This contraint is unfeasable if I limit the Boiler capacity to 2000
subject to Energy_Balance_LT_cstr3 { t in TIME: temp_supply['HPHIGH',t] > 65}:
  Heat_Demand['HPHIGH',t] = sum {c in HEAT_UTIL: c!= 'HPLOW' && c!='HPHIGH'} Heating_HT [c,t];


#Constraints Design Size
subject to SizeTot:
  sum{c in HEAT_UTIL} (Capacity[c]) >= 25000;
subject to SizeMT:
  sum{c in HEAT_UTIL: c!= 'HPLOW'} (Capacity[c])>= 13000;
subject to SizeLT:
  sum{c in HEAT_UTIL: c!= 'HPHIGH'} (Capacity[c])>= 12000;


# Overall Combined heat balance
subject to Energy_Balance_overall_cstr {c in HEAT_UTIL,t in TIME}:
  ComponentSize_t[c,t] = Heating_LT[c,t]+Heating_HT[c,t];

#/!\ DELETED CONSTRAINT

#88% constraint 
subject to EightyeightPerc_Constr:
  sum{h in HP, t in TIME}(ComponentSize_t[h,t]*TIMEsteps[t]) >= 0.88*sum{b in HP,t in TIME}(Heat_Demand[b,t]*TIMEsteps[t]); #SYSTEM REQUIREMENTS
#Fuel Constraint

subject to FUEL_heat_balance_constr{u in FUEL_USERS, t in TIME}: 
 ComponentSize_t[u,t] = FUEL_th_eff[u]*FUEL_Demand[u,t]; #dim #kW

#Elec Balance
subject to FUEL_elec_balance_constr{u in FUEL_USERS, t in TIME}: 
  El_prod[u,t] = FUEL_el_eff[u]*FUEL_Demand[u,t]; #dim kW
subject to El_available_Constr{t in TIME}: #AJOUTER COUTS DES NOUVEAUX PANNEAUX!
  El_prod["SOLAR",t] = ((solarfarm_area+Capacity["SOLAR"])*Efficiency_SolarPanels)*solar_radiation/1000; #dim kW
subject to HP_Energy_Balance_cstr{h in HP,t in TIME}:
  ComponentSize_t[h,t] = COP[h,t]*(-El_prod[h,t]);  #kW

#Total Elec Balance
subject to Elec_demand_system{t in TIME}:
  Elec_Demand[t] = ((spec_annual_elec_demand*floor_area)/(365*24)+ (vol_cooling_water/60 + (Heat_Demand['HPLOW',t]+Heat_Demand['HPHIGH',t])/(4*cp_water*1000))*pumping_cost); #dim #KW 
subject to Electricity_balance_Constr{t in TIME}:
  sum{u in COMPONENTS} El_prod[u,t]+ El_Buy[t] - El_Sell[t]= Elec_Demand[t]; #kW

#Investment
subject to PC_Con{c in COMPONENTS}:
  PC[c] = f_act*((PC_max[c] - PC_min[c])*(Capacity[c]/(C_max[c] - C_min[c])) + Component_Use[c]*PC_min[c]); #USD   # Purchase Cost
# We use a linear interpolation for  the purchase cost between min and max capacity
subject to BM_C_Con{c in COMPONENTS}:
  BM_C[c] = F_P[c]*F_T[c]*F_BM[c]*PC[c];      #Baremoddule Cost  
subject to GR_C_Con{c in COMPONENTS}:
  GR_C[c] = BM_C[c]*(alpha_1*alpha_2 + 1);    #Total cost (Goss cost)
subject to an_CAPEX_Con{c in COMPONENTS}:
  an_CAPEX[c] = GR_C[c]*((interest_rate*(1+interest_rate)^lifetime)/((1+interest_rate)^lifetime - 1));  #Cout annualisés
subject to an_CAPEXTot_Con:
  an_CAPEX_Tot = sum{c in COMPONENTS} an_CAPEX[c];



############################################################################################
# OBJECTIVE FUNCTION
############################################################################################

minimize COST:
sum{t in TIME} ((sum{u in NG_USERS}(c_ng_in*FUEL_Demand[u,t])+ c_ds_in*FUEL_Demand["ICENGINE",t] + c_el_in*El_Buy[t] - c_el_out*El_Sell[t])*TIMEsteps[t])+ an_CAPEX_Tot; #
solve;


############################################################################################
# DISPLAY
############################################################################################

#display COST/(10^6);

display an_CAPEX_Tot/(10^6);
display Heat_Demand;
#display ComponentSize_t;
display Capacity;
#display El_Sell;
#display temp_supply;



end;