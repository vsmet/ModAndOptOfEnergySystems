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
set Base_case;


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
param Year_ind{c in COMPONENTS};
param Ref_ind{c in COMPONENTS};
param F_P{c in COMPONENTS};
param F_T{c in COMPONENTS};
param F_BM {c in COMPONENTS};
param alpha_1;
param alpha_2;
param C_min{c in COMPONENTS};
param C_max{c in COMPONENTS};
param PC_min {c in COMPONENTS};
param PC_max {c in COMPONENTS};  
param f_act{c in COMPONENTS} := Year_ind[c]/Ref_ind[c];
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
param elscen;
param ngscen;

param c_el_in{s in SCENARIO} :=
  if elscen = 1
    then base_case["c_el_in"]*rate[s]
  else base_case["c_el_in"];

param c_el_out{s in SCENARIO} :=
  if elscen = 1
    then base_case["c_el_out"]*rate[s]
  else base_case["c_el_out"];

param c_ng_in{s in SCENARIO} :=
  if ngscen = 1
    then base_case["c_ng_in"]*rate[s]
  else base_case["c_ng_in"];


param c_ds_in{s in SCENARIO} :=
  if ngscen = 1
    then base_case["c_ds_in"]*rate[s]
  else base_case["c_ds_in"];



############################################################################################
# VARIABLES (and defining equations)
############################################################################################

#Variables Investment
var Capacity {c in COMPONENTS, s in SCENARIO} >= 0;
var ElCapacity {c in COMPONENTS, s in SCENARIO} >=0;
var PC{c in COMPONENTS,s in SCENARIO} >= 0;
var BM_C{c in COMPONENTS,s in SCENARIO} >= 0;
var GR_C{c in COMPONENTS,s in SCENARIO} >= 0;
var an_CAPEX {c in COMPONENTS,s in SCENARIO} >=0;
var an_CAPEX_Tot{s in SCENARIO}  >=0;
#var Component_Use {c in COMPONENTS,s in SCENARIO} binary;
var ComponentSize_t {c in COMPONENTS, t in TIME,s in SCENARIO} >= 0;
var ComponentSize_e {c in COMPONENTS, t in TIME,s in SCENARIO} >= 0;
var AnnualCompEnergy {c in COMPONENTS, s in SCENARIO} >= 0;
var AnnualCompEnergyEL {c in COMPONENTS, s in SCENARIO};
var Heating {h in HP,c in HEAT_UTIL, t in TIME,s in SCENARIO} >= 0;
var FUEL_Demand{u in FUEL_USERS, t in TIME,s in SCENARIO}>=0;
#Variables Electricity
var El_prod{u in COMPONENTS, t in TIME,s in SCENARIO};
var El_prod_e{u in COMPONENTS, t in TIME,s in SCENARIO};
var El_Buy{t in TIME,s in SCENARIO} >=0;
var El_Sell{t in TIME,s in SCENARIO}>=0;
var Elec_Demand{t in TIME,s in SCENARIO} >=0;
var Elec_cost{s in SCENARIO};
var Fuel_cost{s in SCENARIO} >= 0;
var Opex{s in SCENARIO};



############################################################################################
# CONSTRAINTS
############################################################################################

#Base case constraints

param Component_Use{c in COMPONENTS,s in SCENARIO}:=
  if c="TURBINE1" || c="TURBINE2" || c="HPHIGH" || c="HPLOW"
    then 1
  else 0;


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


/*
#88% constraint
subject to EightyeightPerc_Constr{s in SCENARIO}:
  sum{h in HP, t in TIME}(ComponentSize_t[h,t,s]*TIMEsteps[t]) >= 0.88*sum{b in HP,t in TIME}(Heat_Demand[b,t]*TIMEsteps[t]); #SYSTEM REQUIREMENTS
*/

#Fuel Constraint
subject to FUEL_heat_balance_constr{u in FUEL_USERS, t in TIME,s in SCENARIO}: 
 ComponentSize_t[u,t,s] = FUEL_th_eff[u]*FUEL_Demand[u,t,s]; 
subject to FUEL_heat_balance_constrCap{u in FUEL_USERS, s in SCENARIO}: 
 (Capacity[u,s]/FUEL_th_eff[u])*FUEL_el_eff[u] = ElCapacity[u,s]; 

#Elec Balance
subject to FUEL_elec_balance_constr{u in FUEL_USERS, t in TIME,s in SCENARIO}: 
  El_prod[u,t,s] = FUEL_el_eff[u]*FUEL_Demand[u,t,s]; #dim kW
subject to El_available_Constr{t in TIME,s in SCENARIO}: #AJOUTER COUTS DES NOUVEAUX PANNEAUX!
  El_prod["SOLAR",t,s] = ((solarfarm_area+Capacity["SOLAR",s])*Efficiency_SolarPanels)*solar_radiation/1000; #dim kW
subject to El_available_Constr2{t in TIME,s in SCENARIO}: #AJOUTER COUTS DES NOUVEAUX PANNEAUX!
  El_prod["SOLAR",t,s] = ElCapacity["SOLAR",s]; #dim kW
subject to HP_Energy_Balance_cstr{h in HP,t in TIME,s in SCENARIO}:
  ComponentSize_t[h,t,s] = COP[h,t]*(-El_prod[h,t,s]);  #kW
subject to HP_Energy_Balance_cstr2{h in HP,s in SCENARIO}:
  Capacity[h,s] = COP[h,5]*(-ElCapacity[h,s]);  #kW

#Total Elec Balance
subject to Elec_demand_system{t in TIME,s in SCENARIO}:
  Elec_Demand[t,s] = ((spec_annual_elec_demand*floor_area)/(365*24)+ (vol_cooling_water/60 + (Heat_Demand['HPLOW',t]+Heat_Demand['HPHIGH',t])/(4*cp_water*1000))*pumping_cost); #dim #KW 
subject to Electricity_balance_Constr{t in TIME,s in SCENARIO}:
  sum{u in COMPONENTS} El_prod[u,t,s]+ El_Buy[t,s] - El_Sell[t,s]= Elec_Demand[t,s]; #kW

#Investment
subject to PC_Con{c in COMPONENTS,s in SCENARIO}:
  PC[c,s] = f_act[c]*((PC_max[c] - PC_min[c])*(Capacity[c,s]/(C_max[c] - C_min[c])) + Component_Use[c,s]*PC_min[c]); #USD   # Purchase Cost
# We use a linear interpolation for  the purchase cost between min and max capacity
subject to BM_C_Con{c in COMPONENTS,s in SCENARIO}:
  BM_C[c,s] = F_P[c]*F_T[c]*F_BM[c]*PC[c,s];      #Baremoddule Cost  
subject to GR_C_Con{c in COMPONENTS,s in SCENARIO}:
  GR_C[c,s] = BM_C[c,s]*(alpha_1*alpha_2 + 1);    #Total cost (Goss cost)
subject to an_CAPEX_Con{c in COMPONENTS,s in SCENARIO}:
  an_CAPEX[c,s] = GR_C[c,s]*((interest_rate*(1+interest_rate)^lifetime)/((1+interest_rate)^lifetime - 1));  #Cout annualisés
subject to an_CAPEXTot_Con{s in SCENARIO}:
  an_CAPEX_Tot[s] = (sum{c in COMPONENTS} an_CAPEX[c,s])/1000;

subject to opex_Con{s in SCENARIO}:
  Opex[s]=sum{t in TIME} ((sum{u in NG_USERS}(c_ng_in[s]*FUEL_Demand[u,t,s])+ c_ds_in[s]*FUEL_Demand["ICENGINE",t,s] + c_el_in[s]*El_Buy[t,s] - c_el_out[s]*El_Sell[t,s])*TIMEsteps[t]);

subject to Power_Energy{c in COMPONENTS, s in SCENARIO, t in TIME}:
  ComponentSize_t[c,t,s]*TIMEsteps[t] = ComponentSize_e[c,t,s];
subject to Energy_Comp{c in COMPONENTS, s in SCENARIO}: #GWh
  (sum{t in TIME}(ComponentSize_e[c,t,s]))/(10^6) = AnnualCompEnergy[c,s];
subject to Fuel_cost_sc{s in SCENARIO}:
  Fuel_cost[s] = (sum{t in TIME} ((sum{u in NG_USERS}(c_ng_in[s]*FUEL_Demand[u,t,s])+ c_ds_in[s]*FUEL_Demand["ICENGINE",t,s] )*TIMEsteps[t]))/1000;
subject to El_cost{s in SCENARIO}:
  Elec_cost[s] = (sum{t in TIME} ((c_el_in[s]*El_Buy[t,s] - c_el_out[s]*El_Sell[t,s])*TIMEsteps[t]))/1000;
subject to Power_EnergyELEC{c in COMPONENTS, s in SCENARIO, t in TIME}:
  El_prod[c,t,s]*TIMEsteps[t] = El_prod_e[c,t,s];
subject to Energy_CompELEC{c in COMPONENTS, s in SCENARIO}: #GWh
  (sum{t in TIME}(El_prod_e[c,t,s]))/(10^6) = AnnualCompEnergyEL[c,s];

#Initialise File
param CapOut, symbolic := "HeatCapacity.csv";
param OpOut, symbolic := "HeatOperation.csv";
param ElOut, symbolic := "ElecBuySell.csv";
param EnOut, symbolic := "HeatEnergy.csv";
param CoOut, symbolic := "Costs.csv";
param ElbOut, symbolic := "ElecCapacity.csv";
param ElEnOut, symbolic := "ElecEnergy.csv";


############################################################################################
# OBJECTIVE FUNCTION
############################################################################################

minimize COST:
sum{s in Base_case}( Opex[s]+1000*an_CAPEX_Tot[s]); #
solve;





############################################################################################
# DISPLAY
############################################################################################
display COST;
display Opex;

printf "Utility_Capacity[kW], " >> CapOut;
for {s in SCENARIO}{
  printf "%s,",s >> CapOut;
}
printf "\n" >> CapOut;
for {c in COMPONENTS: c != "SOLAR"} {
  printf "%s,",c  >> CapOut;
  for {s in SCENARIO}{
    printf "%f,", Capacity[c,s] >> CapOut;
  }
  printf "\n" >> CapOut;
}
printf "end;\n" >> CapOut;

printf "Utility_Capacity[kW], " >> ElbOut;
for {s in SCENARIO}{
  printf "%s,",s >> ElbOut;
}
printf "\n" >> ElbOut;
for {c in COMPONENTS} {
  printf "%s,",c  >> ElbOut;
  for {s in SCENARIO}{
    printf "%f,", ElCapacity[c,s] >> ElbOut;
  }
  printf "\n" >> ElbOut;
}
printf "end;\n" >> ElbOut;

for {s in SCENARIO}{
  printf "%s \n",s >> OpOut;
  # printf "HEATDEMAND," >> OpOut;
  # for {t in TIME}{
  #     printf "%f,",Heat_Demand[t,s]>> OpOut;
  #   }
  # printf "\n" >> OpOut;
  # printf "TIMESTEPS," >> OpOut;
  # for {t in TIME}{
  #     printf "%f,",TIMEsteps[t]>> OpOut;
  #   }
  # printf "\n" >> OpOut;
  printf "\n" >> OpOut;
  for {c in COMPONENTS: c != "SOLAR"}{
    printf "%s,",c >> OpOut;
    for {t in TIME}{
      printf "%f,",ComponentSize_e[c,t,s] >> OpOut;
    }
    printf "\n" >> OpOut;
  }
  printf "\n" >> OpOut;
  printf "\n" >> OpOut;
}




printf "EnergyProduced[GWh], " >> EnOut;
for {s in SCENARIO}{
  printf "%s,",s >> EnOut;
}
printf "\n" >> EnOut;
for {c in COMPONENTS: c!= "SOLAR"} {
  printf "%s,",c  >> EnOut;
  for {s in SCENARIO}{
    printf "%f,", AnnualCompEnergy[c,s] >> EnOut;
  }
  printf "\n" >> EnOut;
}
printf "end;\n" >> EnOut;

printf "EnergyProduced[GWh], " >> ElEnOut;
for {s in SCENARIO}{
  printf "%s,",s >> ElEnOut;
}
printf "\n" >> ElEnOut;
for {c in COMPONENTS} {
  printf "%s,",c  >> ElEnOut;
  for {s in SCENARIO}{
    printf "%f,", AnnualCompEnergyEL[c,s] >> ElEnOut;
  }
  printf "\n" >> ElEnOut;
}
printf "end;\n" >> ElEnOut;




printf "ElectricityBalance [kW] \n" >> ElOut;
for {s in SCENARIO}{
  printf "%s \n",s >> ElOut;
  printf "El_Buy," >> ElOut;
  for {t in TIME}{
    printf "%f,",El_Buy[t,s] >> ElOut;
  }
  printf "\n" >> ElOut;
  printf "El_Sell," >> ElOut;
  for {t in TIME}{
    printf "%f,",El_Sell[t,s] >> ElOut;
  }
  printf "\n" >> ElOut;
  printf "\n" >> ElOut;
}

printf "Costs[kCHF] \n" >> CoOut;
printf "SCENARIO," >> CoOut;
for {s in SCENARIO}{
  printf "%s,",s >> CoOut;
}
printf "\n" >> CoOut;
printf "FUELCOST," >> CoOut;
for {s in SCENARIO}{
  printf "%f,",Fuel_cost[s] >> CoOut;
}
printf "\n" >> CoOut;
printf "ELECCOST," >> CoOut;
for {s in SCENARIO}{
  printf "%f,",Elec_cost[s] >> CoOut;
}
printf "\n" >> CoOut;
printf "CAPEX," >> CoOut;
for {s in SCENARIO}{
  printf "%f,",an_CAPEX_Tot[s] >> CoOut;
}
printf "\n" >> CoOut;

end;