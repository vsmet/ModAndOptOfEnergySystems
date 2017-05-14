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
set SCENARIO;
set BaseC;
set Rate;


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
param interest_rate{s in SCENARIO};
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

#Cost Sensitivity
param base_case{b in BaseC};
param rate{r in Rate};

param c_el_in{s in SCENARIO} :=
  if s="LOW_E" || s="LOW_E_HIGH_ND" || s="LOW_E_LOW_ND"
    then base_case["c_el_in"]*rate["el_low"]
  else if s="HIGH_E"|| s="HIGH_E_HIGH_ND" || s="HIGH_E_LOW_ND"
    then base_case["c_el_in"]*rate["el_high"]
  else base_case["c_el_in"];

param c_el_out{s in SCENARIO} :=
  if s="LOW_E" || s="LOW_E_HIGH_ND" || s="LOW_E_LOW_ND"
    then base_case["c_el_out"]*rate["el_low"]
  else if s="HIGH_E" || s="HIGH_E_HIGH_ND" || s="HIGH_E_LOW_ND"
    then base_case["c_el_out"]*rate["el_high"]
  else base_case["c_el_out"];

param c_ng_in{s in SCENARIO} :=
  if s="LOW_N" || s="LOW_E_LOW_ND" || s="HIGH_E_LOW_ND"
    then base_case["c_ng_in"]*rate["ng_low"]
  else if s="HIGH_N" || s="LOW_E_HIGH_ND" || s="HIGH_E_HIGH_ND"
    then base_case["c_ng_in"]*rate["ng_high"]
  else base_case["c_ng_in"];


param c_ds_in{s in SCENARIO} :=
  if s="LOW_D" || s="LOW_E_LOW_ND" || s="HIGH_E_LOW_ND"
    then base_case["c_ds_in"]*rate["d_low"]
  else if s="HIGH_D" || s="LOW_E_HIGH_ND" || s="HIGH_E_HIGH_ND"
    then base_case["c_ds_in"]*rate["d_high"]
  else base_case["c_ds_in"];



############################################################################################
# VARIABLES (and defining equations)
############################################################################################

#Variables Investment
var Capacity {c in COMPONENTS, s in SCENARIO} >= 0;
var PC{c in COMPONENTS,s in SCENARIO} >= 0;
var BM_C{c in COMPONENTS,s in SCENARIO} >= 0;
var GR_C{c in COMPONENTS,s in SCENARIO} >= 0;
var an_CAPEX {c in COMPONENTS,s in SCENARIO} >=0;
var an_CAPEX_Tot{s in SCENARIO}  >=0;
var Component_Use {c in COMPONENTS,s in SCENARIO} binary;
var ComponentSize_t {c in COMPONENTS, t in TIME,s in SCENARIO} >= 0;
var Heating {h in HP,c in HEAT_UTIL, t in TIME,s in SCENARIO} >= 0;
var FUEL_Demand{u in FUEL_USERS, t in TIME,s in SCENARIO}>=0;
#Variables Electricity
var El_prod{u in COMPONENTS, t in TIME,s in SCENARIO};
var El_Buy{t in TIME,s in SCENARIO} >=0;
var El_Sell{t in TIME,s in SCENARIO}>=0;
var Elec_Demand{t in TIME,s in SCENARIO} >=0;

############################################################################################
# CONSTRAINTS
############################################################################################

#Component Sizing
subject to Size_Constr{c in COMPONENTS, t in TIME,s in SCENARIO}:
  ComponentSize_t[c,t,s] <= Capacity[c,s];             #kW
subject to Component_cmin_cstr {c in COMPONENTS,s in SCENARIO}:
  Component_Use[c,s]*C_min[c] <= Capacity[c,s];
subject to Component_cmax_cstr {c in COMPONENTS,s in SCENARIO}:
  Component_Use[c,s]*C_max[c] >= Capacity[c,s];

#HeatBalance
subject to Heat_dem_Bal{h in HP, t in TIME, s in SCENARIO}:
  Heat_Demand[h,t]= sum{c in HEAT_UTIL}(Heating[h,c,t,s]);
subject to Comp_size_con{c in HEAT_UTIL, t in TIME, s in SCENARIO}:
  ComponentSize_t[c,t,s] =  sum{h in HP}(Heating[h,c,t,s]);
subject to Temp_sup_L{h in HP, t in TIME, s in SCENARIO:temp_supply[h,t] > 50}:
  Heating[h,"HPLOW",t,s] = 0;
subject to Temp_sup_H{t in TIME, s in SCENARIO:temp_supply["HPHIGH",t] > 65}:
  Heating["HPHIGH","HPHIGH",t,s] = 0;

#Constraints Design Size
subject to SizeTot{s in SCENARIO}:
  sum{c in HEAT_UTIL} (Capacity[c,s]) >= 25000;
subject to SizeMT{s in SCENARIO}:
  sum{c in HEAT_UTIL: c!= 'HPLOW'} (Capacity[c,s])>= 13000;
subject to SizeLT{s in SCENARIO}:
  sum{c in HEAT_UTIL: c!= 'HPHIGH'} (Capacity[c,s])>= 12000;


#88% constraint
subject to EightyeightPerc_Constr{s in SCENARIO}:
  sum{h in HP, t in TIME}(ComponentSize_t[h,t,s]*TIMEsteps[t]) >= 0.88*sum{b in HP,t in TIME}(Heat_Demand[b,t]*TIMEsteps[t]); #SYSTEM REQUIREMENTS
#Fuel Constraint
subject to FUEL_heat_balance_constr{u in FUEL_USERS, t in TIME,s in SCENARIO}: 
 ComponentSize_t[u,t,s] = FUEL_th_eff[u]*FUEL_Demand[u,t,s]; #dim #kW

#Elec Balance
subject to FUEL_elec_balance_constr{u in FUEL_USERS, t in TIME,s in SCENARIO}: 
  El_prod[u,t,s] = FUEL_el_eff[u]*FUEL_Demand[u,t,s]; #dim kW
subject to El_available_Constr{t in TIME,s in SCENARIO}: #AJOUTER COUTS DES NOUVEAUX PANNEAUX!
  El_prod["SOLAR",t,s] = ((solarfarm_area+Capacity["SOLAR",s])*Efficiency_SolarPanels)*solar_radiation/1000; #dim kW
subject to HP_Energy_Balance_cstr{h in HP,t in TIME,s in SCENARIO}:
  ComponentSize_t[h,t,s] = COP[h,t]*(-El_prod[h,t,s]);  #kW

#Total Elec Balance
subject to Elec_demand_system{t in TIME,s in SCENARIO}:
  Elec_Demand[t,s] = ((spec_annual_elec_demand*floor_area)/(365*24)+ (vol_cooling_water/60 + (Heat_Demand['HPLOW',t]+Heat_Demand['HPHIGH',t])/(4*cp_water*1000))*pumping_cost); #dim #KW 
subject to Electricity_balance_Constr{t in TIME,s in SCENARIO}:
  sum{u in COMPONENTS} El_prod[u,t,s]+ El_Buy[t,s] - El_Sell[t,s]= Elec_Demand[t,s]; #kW

#Investment
subject to PC_Con{c in COMPONENTS,s in SCENARIO}:
  PC[c,s] = f_act*((PC_max[c] - PC_min[c])*(Capacity[c,s]/(C_max[c] - C_min[c])) + Component_Use[c,s]*PC_min[c]); #USD   # Purchase Cost
# We use a linear interpolation for  the purchase cost between min and max capacity
subject to BM_C_Con{c in COMPONENTS,s in SCENARIO}:
  BM_C[c,s] = F_P[c]*F_T[c]*F_BM[c]*PC[c,s];      #Baremoddule Cost  
subject to GR_C_Con{c in COMPONENTS,s in SCENARIO}:
  GR_C[c,s] = BM_C[c,s]*(alpha_1*alpha_2 + 1);    #Total cost (Goss cost)
subject to an_CAPEX_Con{c in COMPONENTS,s in SCENARIO}:
  an_CAPEX[c,s] = GR_C[c,s]*((interest_rate[s]*(1+interest_rate[s])^lifetime)/((1+interest_rate[s])^lifetime - 1));  #Cout annualisés
subject to an_CAPEXTot_Con{s in SCENARIO}:
  an_CAPEX_Tot[s] = sum{c in COMPONENTS} an_CAPEX[c,s];

#Initialise File
param CapOut, symbolic := "Capacity.csv";


############################################################################################
# OBJECTIVE FUNCTION
############################################################################################

minimize COST:
sum{t in TIME} ((sum{u in NG_USERS}(c_ng_in['HIGH_E_LOW_ND']*FUEL_Demand[u,t,'HIGH_E_LOW_ND'])+ c_ds_in['HIGH_E_LOW_ND']*FUEL_Demand["ICENGINE",t,'HIGH_E_LOW_ND'] + c_el_in['HIGH_E_LOW_ND']*El_Buy[t,'HIGH_E_LOW_ND'] - c_el_out['HIGH_E_LOW_ND']*El_Sell[t,'HIGH_E_LOW_ND'])*TIMEsteps[t])+ an_CAPEX_Tot['HIGH_E_LOW_ND']; #
solve;





############################################################################################
# DISPLAY
############################################################################################
for{c in COMPONENTS}{

  display ComponentSize_t[c,1,"HIGH_E_LOW_ND"];
}
for{h in HP}{
  display temp_supply[h,1];
  display Heat_Demand[h,1];
}
printf "Utility_Capacity[kW], " >> CapOut;
for {s in SCENARIO}{
  printf "%s,",s >> CapOut;
}
printf "\n" >> CapOut;
for {c in COMPONENTS} {
  printf "%s,",c  >> CapOut;
  for {s in SCENARIO}{
    printf "%f,", Capacity[c,s] >> CapOut;
  }
  printf "\n" >> CapOut;
}
printf "end;\n" >> CapOut;


end;