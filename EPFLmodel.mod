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

############################################################################################
# PARAMETER
############################################################################################

/*******************************************************/
# General parameters
/*******************************************************/
param TIMEsteps{t in TIME};					#hr

/*******************************************************/
# Meteo parameters
/*******************************************************/
param external_temp{t in TIME};				#deg C
param solar_radiation{t in TIME};	

/*******************************************************/
# Building parameters
/*******************************************************/
param floor_area{b in BUILDINGS} >= 0;				#m2
param temp_threshold{b in BUILDINGS};				#deg C
param temp_supply_low{b in BUILDINGS,t in TIME} >= 0;	#deg C
param temp_return_low{b in BUILDINGS,t in TIME} >= 0;	#deg C
param temp_supply_high{b in BUILDINGS,t in TIME} >= 0;	#deg C
param temp_return_high{b in BUILDINGS,t in TIME} >= 0;	#deg C

/*******************************************************/
# Demand parameters
/*******************************************************/
param spec_annual_heat_demand{b in BUILDINGS} >= 0, default 0;		#kJ/m2(yr)
param spec_annual_elec_demand{b in BUILDINGS} >= 0, default 0;		#kWh/m2(yr)

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
  Annual_Elec_Demand[b] = floor_area[b] * spec_annual_elec_demand[b];		#kWh(/yr)
  
# Parameter heating signature
param k1{b in BUILDINGS};
param k2{b in BUILDINGS};
 

# TIME-DEPENDENT HEAT DEMAND
var Heat_Demand{b in BUILDINGS, t in TIME} >= 0;
subject to Heat_Demand_Constr{b in BUILDINGS, t in TIME}:
    Heat_Demand[b,t] = 
    if (external_temp[t] < temp_threshold[b]) then
      k1[b] * (external_temp[t]) + k2[b]						#MW
    else 0;
	

# TIME-DEPENDENT ELEC DEMAND
var Elec_Demand{b in BUILDINGS, t in TIME} >= 0;
subject to Elec_Demand_Constr{b in BUILDINGS, t in TIME}:
  Elec_Demand[b,t] = Annual_Elec_Demand[b] / 12;		#kW
 
 
 
 
############################################################################################
# CONSTRAINTS
############################################################################################

# BOILER MODEL
param Efficiency_Boiler := 0.98;

var NG_Demand_Boiler{t in TIME} >= 0;
var Heat_Supple_Boiler{t in TIME} >= 0;
var Capacity_Boiler >= 0;

subject to Boiler_Energy_Balance_Constr{t in TIME}:
  Heat_Supple_Boiler[t] = Efficiency_Boiler*NG_Demand_Boiler[t];	#kW
  
subject to Boiler_Size_Constr{t in TIME}:
  Heat_Supple_Boiler[t] <= Capacity_Boiler;							#kW



#*********** HP MODEL *****************


var EL_Demand_HP1{t in TIME} >=0;
var Heat_Supple_HP1{t in TIME} >= 0;
var Capacity_HP1 >= 0;

subject to HP_Energy_Balance_Constr{t in TIME}:
  Heat_Supple_HP1[t] = COP_low[t]*EL_Demand_HP1[t];  #kW
  
subject to HP1_Size_Constr{t in TIME}:
  Heat_Supple_HP1[t] <= Capacity_HP1;             #kW



var EL_Demand_HP2{t in TIME} >=0;
var Heat_Supple_HP2{t in TIME} >= 0;
var Capacity_HP2 >= 0;

subject to HP_Energy_Balance_Constr{t in TIME}:
  Heat_Supple_HP2[t] = COP_high[t]*EL_Demand_HP2[t];  #kW
  
subject to HP1_Size_Constr{t in TIME}:
  Heat_Supple_HP2[t] <= Capacity_HP2;             #kW  





#SOLAR PANEL MODEL
param Efficiency_SolarPanels := 0.15; #Get a reference !!!!!!!!
param solarfarm_area := 15500; #m^2

var El_Available_Solar{t in TIME} >=0;

subject to El_available_Constr{t in TIME}:
  El_Available_Solar[t] = solarfarm_area*Efficiency_SolarPanels*solar_radiation[t]; #kW
  



# MASS BALANCE NATURAL GAS
var NG_Demand_grid{t in TIME} >= 0;

subject to Natural_gas_Demand_Constr{t in TIME}:
  NG_Demand_grid[t] = NG_Demand_Boiler[t];		#kW
  

# HEAT BALANCE 

subject to EL_balance_Constr{t in TIME}:
  Heat_Supple_HP[t]+Heat_Supple_Boiler[t] = sum{b in BUILDINGS} Heat_Demand[b,t];   #kW


#ELECTRICITY BALANCE

subject to Electricity_balance_Constr{t in TIME}:
  El_Available_Solar[t] + El_Buy[t] - El_Sell[t] - EL_Demand_HP[t]= sum{b in BUILDINGS} Elec_Demand[b,t]; #kW



  
  
  
############################################################################################
# OBJECTIVE FUNCTION
############################################################################################

/*******************************************************/
# Economic parameters
/*******************************************************/
param c_el_in;
param c_el_out;
param c_ng_in;


# To do!
minimize opex:
sum{t in TIME}((c_ng_in*NG_Demand_grid[t] + c_el_in*El_Buy - c_el_out*El_Sell)*TIMEsteps[t]);


solve;

# To do!

display El_Available_Solar;
display El_Buy;
display El_Sell;

end;
