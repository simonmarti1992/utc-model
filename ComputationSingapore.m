%% Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For questions regarding code use, please contact:
%
% Naika Meili (meili@ifu.baug.ethz.ch), 
% ETH Zurich, Future Cities Laboratory, Singapore-ETH Centre, Singapore
%
% Simone Fatichi (fatichi@ifu.baug.ethz.ch), 
% Institute of Environmental Engineering, ETH Zurich, Zurich, Switzerland
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data used in this example was obtained from the authors of:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Velasco E., Roth M., Tan S. H., Quak M., Nabarro S. D. A., and Norford L.(2013).
% The role of vegetation in the CO 2 flux from a tropical urban neighbourhood, 
% Atmospheric Chemistry and Physics, 13, 10 185�10 202, https://doi.org/10.5194/acp-13-10185-2013
%
% Roth M., Jansson C., and Velasco E.(2016). Multi-year energy balance and carbon 
% dioxide fluxes over a residential neighbourhood in a tropical city, 
% International Journal of Climatology, https://doi.org/10.1002/joc.4873
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% RUN TIME SERIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load input parameters
load(fullfile('+data_functions', 'Forcing_Data_Singapore.mat'))

% Assign calculation length
n			=	8760;		% Calculation lenght, e.g. data time series length
m			=	1;			% Usually 1 or length of sensitivity analysis vector

% Assign file names
Name_Site	=	'Singapore';	% Name chosen for "data_functions\Data_UEHM_site_Name.m" 
								% and "data_functions\UEHMForcingData_Name.m"
NameOutput	=	'Singapore_78';	% Name under which the calculation will be saved

% If a new geometrie is chosen, view factors need to be recalculated with
% the MC ray tracing algorithm implemented in UT&C. This will take a couple of
% minutes at the start of the calculation. The results will be saved so further calculations using the
% same urban geometry can just preload the view factors.
OPTION_RAY	=	1; % Load precalculated view factors [1], Recalculate view factors [0]


%% Meteo data: Assign your meteorological input data in belows script MeteoDataRaw using the below specified units(!)
% LWR_in [W/m2]:		Incoming longwave radiation
% SAB1_in [W/m2]:		Incoming direct beam radiation 1
% SAB2_in [W/m2]:		Incoming direct beam radiation 2
% SAD1_in [W/m2]:		Incoming diffuse radiation 1
% SAD2_in [W/m2]:		Incoming diffuse radiation 2
% T_atm	[K]:			Air temperature
% windspeed_u[m/s]:		Wind speed
% pressure_atm [Pa]:	Atmospheric pressure
% rain [mm/h]:			Precipitation
% rel_humidity [-]:		Relative humidity
% Date:					Date and time in Matlabs datetime format.

MeteoDataSG_h.Windspeed(MeteoDataSG_h.Windspeed(1:n,:)==0) = 0.01;	% Wind speed cannot be 0 otherwise the resistance function fails

MeteoDataRaw	=	struct('LWR_in',MeteoDataSG_h.LWRin(1:n,:),'SAB1_in',MeteoDataSG_h.SAB1(1:n,:),...
					'SAB2_in',MeteoDataSG_h.SAB2(1:n,:),'SAD1_in',MeteoDataSG_h.SAD1(1:n,:),...
					'SAD2_in',MeteoDataSG_h.SAD2(1:n,:),'T_atm',MeteoDataSG_h.Tatm(1:n,:),...
					'windspeed_u',MeteoDataSG_h.Windspeed(1:n,:),'pressure_atm',MeteoDataSG_h.Pressure(1:n,:),...
					'rain',MeteoDataSG_h.Precipitation(1:n,:),'rel_humidity',...
					MeteoDataSG_h.RelativeHumidity(1:n,:)./100,'Date',MeteoDataSG_h.date_time(1:n,:));
			


%% Calculation starts here. No need to change anything after this point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
clearvars -except MeteoDataRaw n m Name_Site OPTION_RAY NameOutput

% Meteo data at time step 1 for initialization				
[~,~,HumidityAtm,~,~,~]=feval(strcat('data_functions.UEHMForcingData_',Name_Site),MeteoDataRaw,1);
%[~,~,HumidityAtm,~,~,~]=data_functions.UEHMForcingData(MeteoDataRaw,1);

% Soil parameters
[~,~,~,~,~,~,ParSoilRoof,ParSoilGround,~,~,~,~,~,~,~,~,~,~,~,~]=feval(strcat('data_functions.Data_UEHM_site_',Name_Site),1);
%[~,~,~,~,~,~,ParSoilRoof,ParSoilGround,~,~,~,~,~,~,~,~,~,~,~,~]=data_functions.Data_UEHM_site(1);

ParSoil		=	struct('Roof',ParSoilRoof,'Ground',ParSoilGround);
[~,~,~,~,ParSoil.Roof.O33,~,~,~,~,~]=soil_functions.Soil_parameters(ParSoil.Roof.Psan,ParSoil.Roof.Pcla,ParSoil.Roof.Porg);
[~,~,~,~,ParSoil.Ground.O33,~,~,~,~,~]=soil_functions.Soil_parameters(ParSoil.Ground.Psan,ParSoil.Ground.Pcla,ParSoil.Ground.Porg);
ParSoil.Roof.dz		=	diff(ParSoil.Roof.Zs);	% [mm]  Thickness of the Layers
ParSoil.Ground.dz	=	diff(ParSoil.Ground.Zs);% [mm]  Thickness of the Layers


%% Initializing vectors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Temperature
% TRoofImp		=	Temperature roof impervious area
% TRoofVeg		=	Temperature roof vegetated area
% TRoofIntImp	=	Interior temperature roof impervious area
% TRoofIntVeg	=	Interior temperature roof vegetated area
% TGroundImp	=	Temperature ground impervious area
% TGroundBare	=	Temperature ground bare area
% TGroundVeg	=	Temperature ground vegetated area
% TTree			=	Temperature tree canopy
% TWallSun		=	Temperature sunlit area
% TWallShade	=	Temperature shaded area
% TWallIntSun	=	Interior temperature sunlit wall
% TWallIntShade	=	Interior temperature shaded wall
% TCanyon		=	Temperature canyon
% Tatm			=	Temperature atmosphere(measured)

TempVecNames	=	{'TRoofImp';'TRoofVeg';'TRoofIntImp';'TRoofIntVeg';...
					'TGroundImp';'TGroundBare';'TGroundVeg';'TTree';'TWallSun';...
					'TWallShade';'TWallIntSun';'TWallIntShade';'TCanyon';'Tatm'};

for i=1:size(TempVecNames,1)
	TempVec.(cell2mat(TempVecNames(i)))			=	zeros(n,1,m);
	TempVec.(cell2mat(TempVecNames(i)))(1,:,:)	=	303.16;
end

TempVec.Tatm	=	repmat(MeteoDataRaw.T_atm(:,1),1,1,m); % Temperature atmosphere(measured)

% Dampending temperature
% TGroundImp	=	Dampening temperature ground impervious area
% TGroundBare	=	Dampening temperature ground bare area
% TGroundVeg	=	Dampening temperature ground vegetated area
% TTree			=	Dampening temperature tree canopy

TempDampNames	=	{'TDampGroundImp';'TDampGroundBare';'TDampGroundVeg';'TDampTree'};

for i=1:size(TempDampNames,1)
	TempDamp.(cell2mat(TempDampNames(i)))		=	zeros(n,1,m);
	TempDamp.(cell2mat(TempDampNames(i)))(1,:,:)=	303.16;
end

%% Humidity
HumidityNames	=	{'CanyonRelative';'CanyonSpecific';'CanyonVapourPre';'CanyonRelativeSat';...
					'CanyonSpecificSat';'CanyonVapourPreSat';'AtmRelative';'AtmSpecific';'AtmVapourPre';...
					'AtmRelativeSat';'AtmSpecificSat';'AtmVapourPreSat'};

for i=1:size(HumidityNames,1)
	Humidity.(cell2mat(HumidityNames(i)))	=	zeros(n,1,m);
end
Humidity.CanyonSpecific(1,:,:)				=	HumidityAtm.AtmSpecific;

%% Energy fluxes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Shortwave radiation
% Shortwave radiation absorbed
% SWRabsRoofImp		=	Absorbed shortwave radiation roof impervious area [W/m2 horizontal roof area]
% SWRabsRoofVeg		=	Absorbed shortwave radiation roof vegetated area [W/m2 horizontal roof area]
% SWRabsGroundImp	=	Absorbed shortwave radiation ground impervious area [W/m2 horizontal ground area]
% SWRabsGroundBare	=	Absorbed shortwave radiation ground bare area [W/m2 horizontal ground area]
% SWRabsGroundVeg	=	Absorbed shortwave radiation ground vegetated area [W/m2 horizontal ground area]
% SWRabsTree		=	Absorbed shortwave radiation tree canopy [W/m2 tree sphere??]
% SWRabsWallSun		=	Absorbed shortwave radiation sunlit area [W/m2 vertical wall area]
% SWRabsWallShade	=	Absorbed shortwave radiation shaded area [W/m2 vertical wall area]
% SWRabsTotalRoof	=	Total absorbed shortwave radiation by the roof area [W/m2 horizontal roof area]
% SWRabsTotalGround	=	Total absorbed shortwave radiation by the canyon ground area [W/m2 horizontal ground area]
% SWRabsTotalCanyon	=	Total absorbed shortwave radiation by all the canyon facets [W/m2 horizontal canyon area]
% SWRabsTotalUrban	=	Total absorbed shortwave radiation by all the urban elements (roof plus canyon) [W/m2 horizontal area]
% Incoming shortwave radiation
% SWRinRoofImp		=	Incoming shortwave radiation roof impervious area [W/m2 horizontal roof area]
% SWRinRoofVeg		=	Incoming shortwave radiation roof vegetated area [W/m2 horizontal roof area]
% SWRinGroundImp	=	Incoming shortwave radiation ground impervious area [W/m2 horizontal ground area]
% SWRinGroundBare	=	Incoming shortwave radiation ground bare area [W/m2 horizontal ground area]
% SWRinGroundVeg	=	Incoming shortwave radiation ground vegetated area [W/m2 horizontal ground area]
% SWRinTree			=	Incoming shortwave radiation tree canopy [W/m2 tree sphere??]
% SWRinWallSun		=	Incoming shortwave radiation sunlit area [W/m2 vertical wall area]
% SWRinWallShade	=	Incoming shortwave radiation shaded area [W/m2 vertical wall area]
% SWRinTotalRoof	=	Total incoming shortwave radiation by the roof area [W/m2 horizontal roof area]
% SWRinTotalGround	=	Total incoming shortwave radiation by the canyon ground area [W/m2 horizontal ground area]
% SWRinTotalCanyon	=	Total incoming shortwave radiation by all the canyon facets [W/m2 horizontal canyon area]
% SWRinTotalUrban	=	Total incoming shortwave radiation by all the urban elements (roof plus canyon) [W/m2 horizontal area]
% Outgoing shortwave radiation
% SWRoutRoofImp		=	Outgoing shortwave radiation roof impervious area [W/m2 horizontal roof area]
% SWRoutRoofVeg		=	Outgoing shortwave radiation roof vegetated area [W/m2 horizontal roof area]
% SWRoutGroundImp	=	Outgoing shortwave radiation ground impervious area [W/m2 horizontal ground area]
% SWRoutGroundBare	=	Outgoing shortwave radiation ground bare area [W/m2 horizontal ground area]
% SWRoutGroundVeg	=	Outgoing shortwave radiation ground vegetated area [W/m2 horizontal ground area]
% SWRoutTree		=	Outgoing shortwave radiation tree canopy [W/m2 tree sphere??]
% SWRoutWallSun		=	Outgoing shortwave radiation sunlit area [W/m2 vertical wall area]
% SWRoutWallShade	=	Outgoing shortwave radiation shaded area [W/m2 vertical wall area]
% SWRoutTotalRoof	=	Total outgoing shortwave radiation by the roof area [W/m2 horizontal roof area]
% SWRoutTotalGround	=	Total outgoing shortwave radiation by the canyon ground area [W/m2 horizontal ground area]
% SWRoutTotalCanyon	=	Total outgoing shortwave radiation by all the canyon facets [W/m2 horizontal canyon area]
% SWRoutTotalUrban	=	Total outgoing shortwave radiation by all the urban elements (roof plus canyon) [W/m2 horizontal area]
% Shortwave radiation energy balance
% SWREBRoofImp		=	Energy Balance shortwave radiation roof impervious area [W/m2 horizontal roof area]
% SWREBRoofVeg		=	Energy Balance shortwave radiation roof vegetated area [W/m2 horizontal roof area]
% SWREBGroundImp	=	Energy Balance shortwave radiation ground impervious area [W/m2 horizontal ground area]
% SWREBGroundBare	=	Energy Balance shortwave radiation ground bare area [W/m2 horizontal ground area]
% SWREBGroundVeg	=	Energy Balance shortwave radiation ground vegetated area [W/m2 horizontal ground area]
% SWREBTree			=	Energy Balance shortwave radiation tree canopy [W/m2 tree sphere??]
% SWREBWallSun		=	Energy Balance shortwave radiation sunlit area [W/m2 vertical wall area]
% SWREBWallShade	=	Energy Balance shortwave radiation shaded area [W/m2 vertical wall area]
% SWREBTotalRoof	=	Energy Balance total shortwave radiation by the roof area [W/m2 horizontal roof area]
% SWREBTotalGround	=	Energy Balance total shortwave radiation by the canyon ground area [W/m2 horizontal ground area]
% SWREBTotalCanyon	=	Energy Balance total shortwave radiation by all the canyon facets [W/m2 horizontal canyon area]
% SWREBTotalUrban	=	Energy Balance total outgoing shortwave radiation by all the urban elements (roof plus canyon) [W/m2 horizontal area]

SWRabsNames	=	{'SWRabsRoofImp';'SWRabsRoofVeg';'SWRabsTotalRoof';'SWRabsGroundImp';'SWRabsGroundBare';...
					'SWRabsGroundVeg';'SWRabsTree';'SWRabsWallSun';'SWRabsWallShade';...
					'SWRabsTotalGround';'SWRabsTotalCanyon';'SWRabsTotalUrban'};
				
SWRinNames	=	{'SWRinRoofImp';'SWRinRoofVeg';'SWRinTotalRoof';'SWRinGroundImp';'SWRinGroundBare';...
					'SWRinGroundVeg';'SWRinTree';'SWRinWallSun';'SWRinWallShade';...
					'SWRinTotalGround';'SWRinTotalCanyon';'SWRinTotalUrban'};
								
SWRoutNames	=	{'SWRoutRoofImp';'SWRoutRoofVeg';'SWRoutTotalRoof';'SWRoutGroundImp';'SWRoutGroundBare';...
					'SWRoutGroundVeg';'SWRoutTree';'SWRoutWallSun';'SWRoutWallShade';...
					'SWRoutTotalGround';'SWRoutTotalCanyon';'SWRoutTotalUrban'};
				
SWREBNames	=	{'SWREBRoofImp';'SWREBRoofVeg';'SWREBTotalRoof';'SWREBGroundImp';'SWREBGroundBare';...
					'SWREBGroundVeg';'SWREBTree';'SWREBWallSun';'SWREBWallShade';...
					'SWREBTotalGround';'SWREBTotalCanyon';'SWREBTotalUrban'};

for i=1:size(SWRabsNames,1)
	SWRabs.(cell2mat(SWRabsNames(i)))	=	zeros(n,1,m);
	SWRin.(cell2mat(SWRinNames(i)))		=	zeros(n,1,m);
	SWRout.(cell2mat(SWRoutNames(i)))	=	zeros(n,1,m);
	SWREB.(cell2mat(SWREBNames(i)))		=	zeros(n,1,m);
end

%% Absorbed longwave radiation
% Absorbed longwave radiation
% LWRabsRoofImp		=	Absorbed longwave radiation roof impervious area [W/m2 horizontal roof area]
% LWRabsRoofVeg		=	Absorbed longwave radiation roof vegetated area [W/m2 horizontal roof area]
% LWRabsGroundImp	=	Absorbed longwave radiation ground impervious area [W/m2 horizontal ground area]
% LWRabsGroundBare	=	Absorbed longwave radiation ground bare area [W/m2 horizontal ground area]
% LWRabsGroundVeg	=	Absorbed longwave radiation ground vegetated area [W/m2 horizontal ground area]
% LWRabsTree		=	Absorbed longwave radiation tree canopy [W/m2 tree sphere??]
% LWRabsWallSun		=	Absorbed longwave radiation sunlit area [W/m2 vertical wall area]
% LWRabsWallShade	=	Absorbed longwave radiation shaded area [W/m2 vertical wall area]
% LWRabsTotalRoof	=	Total absorbed longwave radiation by the roof area [W/m2 horizontal roof area]
% LWRabsTotalGround	=	Total absorbed longwave radiation by the canyon ground area [W/m2 horizontal ground area]
% LWRabsTotalCanyon	=	Total absorbed longwave radiation by all the canyon facets [W/m2 horizontal canyon area]
% LWRabsTotalUrban	=	Total absorbed longwave radiation by all the urban elements (roof plus canyon) [W/m2 horizontal area]
% Incoming longwave radiation
% LWRinRoofImp		=	Incoming longwave radiation roof impervious area [W/m2 horizontal roof area]
% LWRinRoofVeg		=	Incoming longwave radiation roof vegetated area [W/m2 horizontal roof area]
% LWRinGroundImp	=	Incoming longwave radiation ground impervious area [W/m2 horizontal ground area]
% LWRinGroundBare	=	Incoming longwave radiation ground bare area [W/m2 horizontal ground area]
% LWRinGroundVeg	=	Incoming longwave radiation ground vegetated area [W/m2 horizontal ground area]
% LWRinTree			=	Incoming longwave radiation tree canopy [W/m2 tree sphere??]
% LWRinWallSun		=	Incoming longwave radiation sunlit area [W/m2 vertical wall area]
% LWRinWallShade	=	Incoming longwave radiation shaded area [W/m2 vertical wall area]
% LWRinTotalRoof	=	Total incoming longwave radiation by the roof area [W/m2 horizontal roof area]
% LWRinTotalGround	=	Total incoming longwave radiation by the canyon ground area [W/m2 horizontal ground area]
% LWRinTotalCanyon	=	Total incoming longwave radiation by all the canyon facets [W/m2 horizontal canyon area]
% LWRinTotalUrban	=	Total incoming longwave radiation by all the urban elements (roof plus canyon) [W/m2 horizontal area]
% Outgoing longwave radiation
% LWRoutRoofImp		=	Outgoing longwave radiation roof impervious area [W/m2 horizontal roof area]
% LWRoutRoofVeg		=	Outgoing longwave radiation roof vegetated area [W/m2 horizontal roof area]
% LWRoutGroundImp	=	Outgoing longwave radiation ground impervious area [W/m2 horizontal ground area]
% LWRoutGroundBare	=	Outgoing longwave radiation ground bare area [W/m2 horizontal ground area]
% LWRoutGroundVeg	=	Outgoing longwave radiation ground vegetated area [W/m2 horizontal ground area]
% LWRoutTree		=	Outgoing longwave radiation tree canopy [W/m2 tree sphere??]
% LWRoutWallSun		=	Outgoing longwave radiation sunlit area [W/m2 vertical wall area]
% LWRoutWallShade	=	Outgoing longwave radiation shaded area [W/m2 vertical wall area]
% LWRoutTotalRoof	=	Total outgoing longwave radiation by the roof area [W/m2 horizontal roof area]
% LWRoutTotalGround	=	Total outgoing longwave radiation by the canyon ground area [W/m2 horizontal ground area]
% LWRoutTotalCanyon	=	Total outgoing longwave radiation by all the canyon facets [W/m2 horizontal canyon area]
% LWRoutTotalUrban	=	Total outgoing longwave radiation by all the urban elements (roof plus canyon) [W/m2 horizontal area]
% Energy Balance of longwave radiation
% LWREBRoofImp		=	Energy Balance longwave radiation roof impervious area [W/m2 horizontal roof area]
% LWREBRoofVeg		=	Energy Balance longwave radiation roof vegetated area [W/m2 horizontal roof area]
% LWREBGroundImp	=	Energy Balance longwave radiation ground impervious area [W/m2 horizontal ground area]
% LWREBGroundBare	=	Energy Balance longwave radiation ground bare area [W/m2 horizontal ground area]
% LWREBGroundVeg	=	Energy Balance longwave radiation ground vegetated area [W/m2 horizontal ground area]
% LWREBTree			=	Energy Balance longwave radiation tree canopy [W/m2 tree sphere??]
% LWREBWallSun		=	Energy Balance longwave radiation sunlit area [W/m2 vertical wall area]
% LWREBWallShade	=	Energy Balance longwave radiation shaded area [W/m2 vertical wall area]
% LWREBTotalRoof	=	Energy Balance total longwave radiation by the roof area [W/m2 horizontal roof area]
% LWREBTotalGround	=	Energy Balance total longwave radiation by the canyon ground area [W/m2 horizontal ground area]
% LWREBTotalCanyon	=	Energy Balance total longwave radiation by all the canyon facets [W/m2 horizontal canyon area]
% LWREBTotalUrban	=	Energy Balance total outgoing longwave radiation by all the urban elements (roof plus canyon) [W/m2 horizontal area]

LWRabsNames	=	{'LWRabsRoofImp';'LWRabsRoofVeg';'LWRabsTotalRoof';'LWRabsGroundImp';'LWRabsGroundBare';...
					'LWRabsGroundVeg';'LWRabsTree';'LWRabsWallSun';'LWRabsWallShade';...
					'LWRabsTotalGround';'LWRabsTotalCanyon';'LWRabsTotalUrban'};
				
LWRinNames	=	{'LWRinRoofImp';'LWRinRoofVeg';'LWRinTotalRoof';'LWRinGroundImp';'LWRinGroundBare';...
					'LWRinGroundVeg';'LWRinTree';'LWRinWallSun';'LWRinWallShade';...
					'LWRinTotalGround';'LWRinTotalCanyon';'LWRinTotalUrban'};
								
LWRoutNames	=	{'LWRoutRoofImp';'LWRoutRoofVeg';'LWRoutTotalRoof';'LWRoutGroundImp';'LWRoutGroundBare';...
					'LWRoutGroundVeg';'LWRoutTree';'LWRoutWallSun';'LWRoutWallShade';...
					'LWRoutTotalGround';'LWRoutTotalCanyon';'LWRoutTotalUrban'};
				
LWREBNames	=	{'LWREBRoofImp';'LWREBRoofVeg';'LWREBTotalRoof';'LWREBGroundImp';'LWREBGroundBare';...
					'LWREBGroundVeg';'LWREBTree';'LWREBWallSun';'LWREBWallShade';...
					'LWREBTotalGround';'LWREBTotalCanyon';'LWREBTotalUrban'};

for i=1:size(LWRabsNames,1)
	LWRabs.(cell2mat(LWRabsNames(i)))	=	zeros(n,1,m);
	LWRin.(cell2mat(LWRinNames(i)))		=	zeros(n,1,m);
	LWRout.(cell2mat(LWRoutNames(i)))	=	zeros(n,1,m);
	LWREB.(cell2mat(LWREBNames(i)))		=	zeros(n,1,m);
end

%% Sensible heat flux
% HfluxRoofImp		=	Sensible heat flux of impervious roof area to atmosphere [W/m2 horizontal roof area]
% HfluxRoofVeg		=	Sensible heat flux of vegetated roof area to atmosphere [W/m2 horizontal roof area]
% HfluxGroundImp	=	Sensible heat flux of impervious ground area to canyon [W/m2 horizontal ground area]
% HfluxGroundBare	=	Sensible heat flux of bare ground area to canyon [W/m2 horizontal ground area]
% HfluxGroundVeg	=	Sensible heat flux of vegetated ground area to canyon [W/m2 horizontal ground area]
% HfluxGround		=	Sensible heat flux of impground area to canyon [W/m2 horizontal ground area]
% HfluxTree			=	Sensible heat flux of tree canopy to canyon [W/m2 tree sphere??]
% HfluxWallSun		=	Sensible heat flux of sunlit wall to canyon [W/m2 vertical wall area]
% HfluxWallShade	=	Sensible heat flux of shaded wall to canyon [W/m2 vertical wall area]
% HfluxCanyon		=	Sensible heat flux of canyon to atmosphere [W/m2 horizontal canyon area]
% HfluxRoof			=	Total sensible heat flux of roof area to atmosphere [W/m2 horizontal roof area]
% HfluxUrban		=	Total sensible heat flux of urban area to atmosphere [W/m2 horizontal area]
				
HfluxNames	=	{'HfluxRoofImp';'HfluxRoofVeg';'HfluxRoof';'HfluxGroundImp';'HfluxGroundBare';...
					'HfluxGroundVeg';'HfluxGround';'HfluxTree';'HfluxWallSun';'HfluxWallShade';...
					'HfluxCanyon';'HfluxUrban'};

for i=1:size(HfluxNames,1)
	Hflux.(cell2mat(HfluxNames(i)))	=	zeros(n,1,m);
end

%% Latent heat flux
% LEfluxRoofImp			=	Latent heat flux of intercepted water from impervious roof area to atmosphere (LEroof_imp_pond) [W/m2 horizontal roof area]
% LEfluxRoofVegInt		=	Latent heat flux of intercepted water on roof vegetation to atmosphere (LEroof_veg_int) [W/m2 horizontal roof area]
% LEfluxRoofVegPond		=	Latent heat flux of intercepted water on ground under roof vegetation to atmosphere (LEroof_veg_pond) [W/m2 horizontal roof area]
% LEfluxRoofVegSoil		=	Latent heat flux of water from roof soil under vegetation to atmosphere (LEroof_veg_soil) [W/m2 horizontal roof area]
% LTEfluxRoofVeg		=	Latent heat flux of transpiration from roof plants to atmosphere (LTEroof_veg) [W/m2 horizontal roof area]
% LEfluxRoofVeg			=	Total latent heat flux of vegetated roof to atmosphere [W/m2 horizontal roof area]
% LEfluxRoof			=	Total latent heat flux of roof to atmosphere [W/m2 horizontal roof area]
% LEfluxGroundImp		=	Latent heat flux of intercepted water on impervious ground area to canyon (LEground_imp_pond)[W/m2 horizontal ground area]
% LEfluxGroundBarePond	=	Latent heat flux of  water on bare ground to canyon (LEground_bare_pond)[W/m2 horizontal ground area]
% LEfluxGroundBareSoil	=	Latent heat flux of  water from bare ground to canyon (LEground_bare_soil) [W/m2 horizontal ground area]
% LEfluxGroundBare		=	Total latent heat flux of bare ground area to canyon [W/m2 horizontal ground area]
% LEfluxGroundVegInt	=	Latent heat flux of intercepted water on ground vegetation to canyon (LEground_veg_int) [W/m2 horizontal ground area]
% LEfluxGroundVegPond	=	Latent heat flux of intercepted water on ground under vegetation to canyon (LEground_veg_pond) [W/m2 horizontal ground area]
% LEfluxGroundVegSoil	=	Latent heat flux of water from ground soil under vegetation to canyon (LEground_veg_soil) [W/m2 horizontal ground area]
% LTEfluxGroundVeg		=	Latent heat flux of transpiration from ground plants to canyon (LTEground_veg) [W/m2 horizontal ground area]
% LEfluxGroundVeg		=	Total latent heat flux of vegetated ground to canyon [W/m2 horizontal ground area]
% LEfluxGround			=	Total latent heat flux of ground to canyon [W/m2 horizontal roof area]
% LEfluxTreeInt			=	Latent heat flux of intercepted water on tree canopy to canyon (LE_tree_int) [W/m2 tree sphere??]
% LTEfluxTree			=	Latent heat flux of transpiration from tree canopy to canyon (LTE_tree) [W/m2 tree sphere??]
% LEfluxTree			=	Total latent heat flux of tree canopy to canyon [W/m2 tree sphere??]
% LEfluxWallSun			=	Latent heat flux of sunlit wall to canyon [W/m2 vertical wall area]
% LEfluxWallShade		=	Latent heat flux of shaded wall to canyon [W/m2 vertical wall area]
% LEfluxCanyon			=	Latent heat flux of canyon to atmosphere [W/m2 horizontal canyon area]
% LEfluxUrban			=	Total latent heat flux of urban area to atmosphere [W/m2 horizontal area]

LEfluxNames	=	{'LEfluxRoofImp';'LEfluxRoofVegInt';'LEfluxRoofVegPond';'LEfluxRoofVegSoil';...
					'LTEfluxRoofVeg';'LEfluxRoofVeg';'LEfluxRoof';'LEfluxGroundImp';'LEfluxGroundBarePond';...
					'LEfluxGroundBareSoil';'LEfluxGroundBare';'LEfluxGroundVegInt';...
					'LEfluxGroundVegPond';'LEfluxGroundVegSoil';'LTEfluxGroundVeg';'LEfluxGroundVeg';...
					'LEfluxGround';'LEfluxTreeInt';'LTEfluxTree';'LEfluxTree';'LEfluxWallSun';...
					'LEfluxWallShade';'LEfluxCanyon';'LEfluxUrban'};

for i=1:size(LEfluxNames,1)
	LEflux.(cell2mat(LEfluxNames(i)))	=	zeros(n,1,m);
end

%% Conductive heat fluxes
% G1RoofImp		=	Conductive heat flux of first layer of impervious roof
% G1RoofVeg		=	Conductive heat flux of first layer of vegetated roof
% G2RoofImp		=	Conductive heat flux of second layer of impervious roof
% G2RoofVeg		=	Conductive heat flux of second layer of vegetated roof
% G1GroundImp	=	Conductive heat flux of impervious ground (G_groundimp)
% G1GroundBare	=	Conductive heat flux of bare ground (G_groundbare)
% G1GroundVeg	=	Conductive heat flux of vegetated ground (G_groundveg)
% GTree			=	Conductive heat flux tree
% G1WallSun		=	Conductive heat flux of first layer of sunlit wall (G1_wallsun)
% G1WallShade	=	Conductive heat flux of first layer of shaded wall (G1_wallshade)
% G2WallSun		=	Conductive heat flux of second layer of sunlit wall (G2_wallsun)
% G2WallShade	=	Conductive heat flux of second layer of shaded wall (G2_wallshade)
	
GfluxNames	=	{'G1RoofImp';'G1RoofVeg';'G2RoofImp';'G2RoofVeg';'G1Roof';'G2Roof';...
					'G1GroundImp';'G1GroundBare';'G1GroundVeg';'G1Ground';'GTree';'G1WallSun';...
					'G1WallShade';'G2WallSun';'G2WallShade';'G1Canyon';'G2Canyon';'G1Urban';'G2Urban'};

for i=1:size(GfluxNames,1)
	Gflux.(cell2mat(GfluxNames(i)))	=	zeros(n,1,m);
end

%% Heat storage in surfaces
% dsRoofImp		=	Storage of energy in impervious roof
% dsRoofVeg		=	Storage of energy in vegetated roof
% dsGroundImp	=	Storage of energy in impervious ground
% dsGroundBare	=	Storage of energy in bare ground
% dsGroundVeg	=	Storage of energy in vegetated ground
% dsTree		=	Storage of energy in tree canopy
% dsWallSun		=	Storage of energy in sunlit wall 
% dsWallShade	=	Storage of energy in shaded wall
% dsCanyonAir	=	Storage of energy in canyon air

dStorageNames	=	{'dsRoofImp';'dsRoofVeg';'dsRoof';'dsGroundImp';'dsGroundBare';...
					'dsGroundVeg';'dsTree';'dsWallSun';'dsWallShade';'dsCanyonAir'};

for i=1:size(dStorageNames,1)
	dStorage.(cell2mat(dStorageNames(i)))	=	zeros(n,1,m);
end

%% Resistances
% raRooftoAtm	=	Aerodynamic resistance ra from roof to atmosphere
% rap_LRoof		=	Undercanopy resistance rap_L roof
% rb_LRoof		=	Leaf boundary resistance rb_L roof
% r_soilRoof	=	Soil resistance rb_soil roof
% rs_sunRoof	=	Stomata resistance sunlit vegetation rs_sun_roof
% rs_shdRoof	=	Stomata resistance shaded vegetation rs_shd_roof
% raCanyontoAtm	=	Aerodynamic resistance ra from canyon to atmosphere
% raGroundtoAtm	=	Aerodynamic resistance ra from ground to canyon air
% raTreetoAtm	=	Aerodynamic resistance ra from tree to canyon air
% raWalltoAtm	=	Aerodynamic resistance ra from wall to canyon air
% rap_HGround	=	Undercanopy resistance rap_H ground
% rap_LGround	=	Undercanopy resistance rap_L ground
% rb_HGround	=	Leaf boundary resistance rb_H ground
% rb_LGround	=	Leaf boundary resistance rb_L ground
% r_soilGround	=	Soil resistance rb_soil ground
% rs_sunGround	=	Stomata resistance sunlit vegetation rs_sun_ground
% rs_shdGround	=	Stomata resistance shaded vegetation rs_shd_ground
% rs_sunTree	=	Stomata resistance sunlit vegetation rs_sun_tree
% rs_shdTree	=	Stomata resistance shaded vegetation rs_shd_ground
	
RESNames	=	{'raRooftoAtm';'rap_LRoof';'rb_LRoof';'r_soilRoof';'rs_sunRoof';'rs_shdRoof';...
					'raCanyontoAtm';'rap_can';'rap_Htree_In';'rb_HGround';'rb_LGround';...
					'r_soilGroundbare';'r_soilGroundveg';'alp_soilGroundbare';'alp_soilGroundveg';...
					'rs_sunGround';'rs_shdGround';'rs_sunTree';'rs_shdTree';...
					'RES_w1';'RES_w2';'rap_W1_In';'rap_W2_In'};

for i=1:size(RESNames,1)
	RES.(cell2mat(RESNames(i)))	=	zeros(n,1,m);
end

%% Water fluxes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Evapotranspiration
% EfluxRoofImp			=	Evaporation flux of intercepted water from impervious roof area to atmosphere (Eroof_imp_pond) [W/m2 horizontal roof area]
% EfluxRoofVegInt		=	Evaporation flux of intercepted water on roof vegetation to atmosphere (Eroof_veg_int) [W/m2 horizontal roof area]
% EfluxRoofVegPond		=	Evaporation flux of intercepted water on ground under roof vegetation to atmosphere (Eroof_veg_pond) [W/m2 horizontal roof area]
% EfluxRoofVegSoil		=	Evaporation flux of water from roof soil under vegetation to atmosphere (Eroof_veg_soil) [W/m2 horizontal roof area]
% TEfluxRoofVeg			=	Evaporation flux of transpiration from roof plants to atmosphere (TEroof_veg) [W/m2 horizontal roof area]
% EfluxRoofVeg			=	Total evaporation flux of vegetated roof to atmosphere [W/m2 horizontal roof area]
% EfluxRoof				=	Total evaporation flux of roof to atmosphere [W/m2 horizontal roof area]
% EfluxGroundImp		=	Evaporation flux of intercepted water on impervious ground area to canyon (Eground_imp_pond)[W/m2 horizontal ground area]
% EfluxGroundBarePond	=	Evaporation flux of  water on bare ground to canyon (Eground_bare_pond)[W/m2 horizontal ground area]
% EfluxGroundBareSoil	=	Evaporation flux of  water from bare ground to canyon (Eground_bare_soil) [W/m2 horizontal ground area]
% EfluxGroundBare		=	Total evaporation flux of bare ground area to canyon [W/m2 horizontal ground area]
% EfluxGroundVegInt		=	Evaporation flux of intercepted water on ground vegetation to canyon (Eground_veg_int) [W/m2 horizontal ground area]
% EfluxGroundVegPond	=	Evaporation flux of intercepted water on ground under vegetation to canyon (Eground_veg_pond) [W/m2 horizontal ground area]
% EfluxGroundVegSoil	=	Evaporation flux of water from ground soil under vegetation to canyon (Eground_veg_soil) [W/m2 horizontal ground area]
% TEfluxGroundVeg		=	Evaporation flux of transpiration from ground plants to canyon (TEground_veg) [W/m2 horizontal ground area]
% EfluxGroundVeg		=	Total evaporation flux of vegetated ground to canyon [W/m2 horizontal ground area]
% EfluxGround			=	Total evaporation flux of ground to canyon [W/m2 horizontal roof area]
% EfluxTreeInt			=	Evaporation flux of intercepted water on tree canopy to canyon (E_tree_int) [W/m2 tree sphere??]
% TEfluxTree			=	Evaporation flux of transpiration from tree canopy to canyon (TE_tree) [W/m2 tree sphere??]
% EfluxTree				=	Total evaporation flux of tree canopy to canyon [W/m2 tree sphere??]
% EfluxWallSun			=	Evaporation flux of sunlit wall to canyon [W/m2 vertical wall area]
% EfluxWallShade		=	Evaporation flux of shaded wall to canyon [W/m2 vertical wall area]
% EfluxCanyon			=	Evaporation flux of canyon to atmosphere [W/m2 horizontal canyon area]
% EfluxUrban			=	Total evaporation flux of urban area to atmosphere [W/m2 horizontal area]

EfluxNames	=	{'EfluxRoofImp';'EfluxRoofVegInt';'EfluxRoofVegPond';'EfluxRoofVegSoil';...
					'TEfluxRoofVeg';'EfluxRoofVeg';'EfluxRoof';'EfluxGroundImp';'EfluxGroundBarePond';...
					'EfluxGroundBareSoil';'EfluxGroundBare';'EfluxGroundVegInt';...
					'EfluxGroundVegPond';'EfluxGroundVegSoil';'TEfluxGroundVeg';'EfluxGroundVeg';...
					'EfluxGround';'EfluxTreeInt';'TEfluxTree';'EfluxTree';'EfluxWallSun';...
					'EfluxWallShade';'EfluxCanyon';'EfluxUrban'};

for i=1:size(EfluxNames,1)
	Eflux.(cell2mat(EfluxNames(i)))	=	zeros(n,1,m);
end

%% Runoff / Runon
% QRoofImp			=	Runoff of impervious area of roof (q_runon_imp)
% QRoofVegDrip		=	Runoff, Dripping, etc from ground vegetation to roof ground (q_runon_veg)
% QRoofVegPond		=	Runoff from roof ground under vegetation due to limitation in infiltration capacity (q_runon_ground_veg)
% QRoofVegSoil		=	Runoff due to roof soil saturation (Rd_veg)
% QGroundImp		=	Runoff of impervious area of ground (q_runon_imp)
% QGroundBarePond	=	Runoff of bare area of ground due to limitation in infiltration capacity(q_runon_bare)
% QGroundBareSoil	=	Runoff of bare area of ground due to soil saturation (Rd_bare)
% QTree				=	Runoff, Dripping, etc from tree to ground (q_runon_tree)
% QGroundVegDrip	=	Runoff, Dripping, etc from ground vegetation to roof ground (q_runon_veg)
% QGroundVegPond	=	Runoff from roof ground under vegetation due to limitation in infiltration capacity (q_runon_ground_veg)
% QGroundVegSoil	=	Runoff due to roof soil saturation (Rd_veg)

RunoffNames	=	{'QRoofImp';'QRoofVegDrip';'QRoofVegPond';'QRoofVegSoil';...
					'QGroundImp';'QGroundBarePond';'QGroundBareSoil';'QTree';'QGroundVegDrip';...
					'QGroundVegPond';'QGroundVegSoil'};

for i=1:size(RunoffNames,1)
	Runoff.(cell2mat(RunoffNames(i)))	=	zeros(n,1,m);
end

% RunonRoofTot		=	Total roof runon to the next time step
% RunoffRoofTot		=	Total roof runoff that is removed from the system
% RunonGroundTot	=	Total ground runon to the next time step
% RunoffGroundTot	=	Total ground runoff that is removed from the system

RunonNames	=	{'RunonRoofTot';'RunoffRoofTot';'RunonGroundTot';'RunoffGroundTot';'RunonUrban';'RunoffUrban'};

for i=1:size(RunonNames,1)
	Runon.(cell2mat(RunonNames(i)))	=	zeros(n,1,m);
end

%% Leakage
% LkRoofImp		=	Leakage from impervious roof (Lk_imp)
% LkRoofVeg		=	Leakage from last soil layer of vegetated roof (Lk_soil_veg)
% LkRoof			=	Total leakage of roof
% LkGroundImp		=	Leakage from impervious ground (Lk_imp)
% LkGroundBare	=	Leakage from last soil layer of bare ground (Lk_soil_bare)
% LkGroundVeg		=	Leakage from last soil layer of vegetated ground (Lk_soil_veg)
% LkGround		=	Total leakage of ground

LeakageNames	=	{'LkRoofImp';'LkRoofVeg';'LkRoof';'LkGroundImp';...
					'LkGroundBare';'LkGroundVeg';'LkGround';'LkUrban'};

for i=1:size(LeakageNames,1)
	Leakage.(cell2mat(LeakageNames(i)))	=	zeros(n,1,m);
end

%% Interception / change in interception
% IntRoofImp		=	Interception on impervious roof area (In_ground_imp)
% IntRoofVegPlant	=	Interception on plant surfaces (In_ground_veg)
% IntRoofVegGround	=	Interception on ground (In_ground_underveg)
% IntGroundImp		=	Interception on impervious ground area (In_ground_imp)
% IntGroundBare		=	Interception on bare ground area (In_ground_bare)
% IntGroundVegPlant	=	Interception on plant surfaces (In_ground_veg)
% IntGroundVegGround=	Interception on ground (In_ground_underveg)
% IntTree			=	Interception on tree (In_tree)

IntNames	=	{'IntRoofImp';'IntRoofVegPlant';'IntRoofVegGround';'IntRooftot';'IntGroundImp';...
					'IntGroundBare';'IntGroundVegPlant';'IntGroundVegGround';'IntTree'};

for i=1:size(IntNames,1)
	Int.(cell2mat(IntNames(i)))	=	zeros(n,1,m);
end

% Change in interception
% dInt_dtRoofImp		=	Change in interception on impervious roof area (dIn_imp_dt)
% dInt_dtRoofVegPlant	=	Change in interception on plant surfaces (dIn_veg_dt)
% dInt_dtRoofVegGround	=	Change in interception on ground (dIn_ground_veg_dt)
% dInt_dtGroundImp		=	Change in interception on impervious ground area (dIn_imp_dt)
% dInt_dtGroundBare		=	Change in interception on bare ground area (dIn_bare_dt)
% dInt_dtGroundVegPlant	=	Change in interception on plant surfaces (dIn_veg_dt)
% dInt_dtGroundVegGround=	Change in interception on ground (dIn_ground_veg_dt)
% dInt_dtTree			=	Change in interception on tree (dIn_tree_dt)

dInt_dtNames	=	{'dInt_dtRoofImp';'dInt_dtRoofVegPlant';'dInt_dtRoofVegGround';'dInt_dtRooftot';'dInt_dtGroundImp';...
					'dInt_dtGroundBare';'dInt_dtGroundVegPlant';'dInt_dtGroundVegGround';'dInt_dtTree'};

for i=1:size(dInt_dtNames,1)
	dInt_dt.(cell2mat(dInt_dtNames(i)))	=	zeros(n,1,m);
end

%% Infiltration
% fRoofVeg	=	Infiltration in first soil layer of vegetated roof (f_roof_veg)
% fGroundBare	=	Infiltration in first soil layer of bare ground (f_ground_bare)
% fGroundVeg	=	Infiltration in first soil layer of vegetated ground (f_ground_veg)	

InfiltrationNames	=	{'fRoofVeg';'fGroundBare';'fGroundVeg';'fGroundImp'};

for i=1:size(InfiltrationNames,1)
	Infiltration.(cell2mat(InfiltrationNames(i)))	=	zeros(n,1,m);
end

%% Water volume in soil / change in water volume in soil
% Initializing soil water content in the first time step.
% I chose field capacity O33 as a starting point
Vwater						=	[];
Vwater.VRoofSoilVeg			=	zeros(n,ParSoil.Roof.ms,m);		%  Water volume in the different soil layers of roof (Vw_soil)
Vwater.VGroundSoilImp		=	zeros(n,ParSoil.Ground.ms,m);	%  Water volume in the different soil layers of ground (Vw_soil)
Vwater.VGroundSoilBare		=	zeros(n,ParSoil.Ground.ms,m);	%  Water volume in the different soil layers of ground (Vw_soil)
Vwater.VGroundSoilVeg		=	zeros(n,ParSoil.Ground.ms,m);	%  Water volume in the different soil layers of ground (Vw_soil)
Vwater.VGroundSoilTot		=	zeros(n,ParSoil.Ground.ms,m);	%  Water volume in the different soil layers of ground (Vw_soil)

Vwater.VRoofSoilVeg(1,:,:)	=	repmat(ParSoil.Roof.O33.*ParSoil.Roof.dz,1,1,m);		% Starting point at field capacity
Vwater.VGroundSoilImp(1,:,:)=	repmat(ParSoil.Ground.O33.*ParSoil.Ground.dz,1,1,m);	% Starting point at field capacity
Vwater.VGroundSoilBare(1,:,:)=	repmat(ParSoil.Ground.O33.*ParSoil.Ground.dz,1,1,m);	% Starting point at field capacity
Vwater.VGroundSoilVeg(1,:,:)=	repmat(ParSoil.Ground.O33.*ParSoil.Ground.dz,1,1,m);	% Starting point at field capacity
Vwater.VGroundSoilTot(1,:,:)=	repmat(ParSoil.Ground.O33.*ParSoil.Ground.dz,1,1,m);	% Starting point at field capacity

dVwater_dt						=	[];
dVwater_dt.dVRoofSoilVeg_dt		=	zeros(n,1,m);	% Maybe only total? I should check zeros(n,1);
dVwater_dt.dVGroundSoilImp_dt	=	zeros(n,1,m); % Maybe only total? I should check zeros(n,1);
dVwater_dt.dVGroundSoilBare_dt	=	zeros(n,1,m); % Maybe only total? I should check zeros(n,1);
dVwater_dt.dVGroundSoilVeg_dt	=	zeros(n,1,m); % Maybe only total? I should check zeros(n,1);
dVwater_dt.dVGroundSoilTot_dt	=	zeros(n,1,m); % Maybe only total? I should check zeros(n,1);

% dVwater_dt						=	[];
% dVwater_dt.dVRoofSoilVeg_dt		=	zeros(n,ParSoil.Roof.ms,m);	% Maybe only total? I should check zeros(n,1);
% dVwater_dt.dVGroundSoilImp_dt	=	zeros(n,ParSoil.Ground.ms,m); % Maybe only total? I should check zeros(n,1);
% dVwater_dt.dVGroundSoilBare_dt	=	zeros(n,ParSoil.Ground.ms,m); % Maybe only total? I should check zeros(n,1);
% dVwater_dt.dVGroundSoilVeg_dt	=	zeros(n,ParSoil.Ground.ms,m); % Maybe only total? I should check zeros(n,1);
% dVwater_dt.dVGroundSoilTot_dt	=	zeros(n,ParSoil.Ground.ms,m); % Maybe only total? I should check zeros(n,1);

%% Soil moisture
Owater						=	[];
Owater.OwRoofSoilVeg		=	zeros(n,ParSoil.Roof.ms,m);		%  Soil moisture in the different soil layers of roof
Owater.OwGroundSoilImp		=	zeros(n,ParSoil.Ground.ms,m);		%  Soil moisture in the different soil layers of ground
Owater.OwGroundSoilBare		=	zeros(n,ParSoil.Ground.ms,m);		%  Soil moisture in the different soil layers of ground
Owater.OwGroundSoilVeg		=	zeros(n,ParSoil.Ground.ms,m);		%  Soil moisture in the different soil layers of ground
Owater.OwGroundSoilTot		=	zeros(n,ParSoil.Ground.ms,m);		%  Soil moisture in the different soil layers of ground

Owater.OwRoofSoilVeg(1,:,:)	=	ParSoil.Roof.O33;				% Starting point at field capacity
Owater.OwGroundSoilImp(1,:,:)=	ParSoil.Ground.O33;				% Starting point at field capacity
Owater.OwGroundSoilBare(1,:,:)=	ParSoil.Ground.O33;				% Starting point at field capacity
Owater.OwGroundSoilVeg(1,:,:)=	ParSoil.Ground.O33;				% Starting point at field capacity
Owater.OwGroundSoilTot(1,:,:)=	ParSoil.Ground.O33;				% Starting point at field capacity

Owater.OwGroundSoilImp(:,1:2,:)=	NaN;				% Starting point at field capacity

OSwater						=	[];
OSwater.OSwRoofSoilVeg		=	zeros(n,ParSoil.Roof.ms,m);
OSwater.OSwGroundSoilImp	=	zeros(n,ParSoil.Ground.ms,m);
OSwater.OSwGroundSoilBare	=	zeros(n,ParSoil.Ground.ms,m);
OSwater.OSwGroundSoilVeg	=	zeros(n,ParSoil.Ground.ms,m);
OSwater.OSwGroundSoilTot	=	zeros(n,ParSoil.Ground.ms,m);

% Lateral soil water flux
Qinlat				=	[];
Qinlat.Qin_bare2imp	=	zeros(n,ParSoil.Ground.ms,m);
Qinlat.Qin_veg2imp	=	zeros(n,ParSoil.Ground.ms,m);
Qinlat.Qin_veg2bare	=	zeros(n,ParSoil.Ground.ms,m);
Qinlat.Qin_imp2bare	=	zeros(n,ParSoil.Ground.ms,m);
Qinlat.Qin_bare2veg	=	zeros(n,ParSoil.Ground.ms,m);
Qinlat.Qin_imp2veg	=	zeros(n,ParSoil.Ground.ms,m);
Qinlat.Qin_imp		=	zeros(n,ParSoil.Ground.ms,m);
Qinlat.Qin_bare		=	zeros(n,ParSoil.Ground.ms,m);
Qinlat.Qin_veg		=	zeros(n,ParSoil.Ground.ms,m);

%% Max extractable water and soil water potential for plants in soil
% Max extractable water
% Extractable water: Soil moisture in the different soil layers (Exwat)
ExWaterNames	=	{'ExWaterRoofVeg_H';'ExWaterRoofVeg_L';...
					'ExWaterGroundImp_H';'ExWaterGroundImp_L';...
					'ExWaterGroundBare_H';'ExWaterGroundBare_L';...
					'ExWaterGroundVeg_H';'ExWaterGroundVeg_L';...
					'ExWaterGroundTot_H';'ExWaterGroundTot_L'};
for i=1:2
	ExWater.(cell2mat(ExWaterNames(i)))	=	zeros(n,ParSoil.Roof.ms,m);
end

for i=3:size(ExWaterNames,1)
	ExWater.(cell2mat(ExWaterNames(i)))	=	zeros(n,ParSoil.Ground.ms,m);
end

% Soil water potential for plants in soil
% SoilPotWRoof_H	=	soil water potential for plants (Psi_s_H)
% SoilPotWRoof_L	=	soil water potential for plants (Psi_s_L)
% SoilPotWGround_H	=	soil water potential for plants (Psi_s_H)
% SoilPotWGround_L	=	soil water potential for plants (Psi_s_L)
SoilPotWNames	=	{'SoilPotWRoofVeg_H';'SoilPotWRoofVeg_L';...
					'SoilPotWGroundImp_H';'SoilPotWGroundImp_L';...
					'SoilPotWGroundBare_H';'SoilPotWGroundBare_L';...
					'SoilPotWGroundVeg_H';'SoilPotWGroundVeg_L';...
					'SoilPotWGroundTot_H';'SoilPotWGroundTot_L'};

for i=1:size(SoilPotWNames,1)
	SoilPotW.(cell2mat(SoilPotWNames(i)))	=	zeros(n,1,m);
end

%% Ci = Leaf Interior  CO2 concentration [umolCO2/mol]
% CiCO2LeafGroundVegSun	=	Ci_sun_veg
% CiCO2LeafGroundVegShd	=	Ci_shd_veg
% CiCO2LeafTreeSun		=	Ci_sun_tree
% CiCO2LeafTreeShd		=	Ci_shd_tree

CiCO2LeafNames	=	{'CiCO2LeafRoofVegSun';'CiCO2LeafRoofVegShd';...
					'CiCO2LeafGroundVegSun';'CiCO2LeafGroundVegShd';...
					'CiCO2LeafTreeSun';'CiCO2LeafTreeShd'};

for i=1:size(CiCO2LeafNames,1)
	CiCO2Leaf.(cell2mat(CiCO2LeafNames(i)))			=	zeros(n,1,m);
	CiCO2Leaf.(cell2mat(CiCO2LeafNames(i)))(1,:,:)	=	400;
end

%% Energy and Water Balance
WBRoofNames	=	{'WBRoofImp';'WBRoofVegInVeg';'WBRoofVegInGround';'WBRoofVegSoil';...
				'WBRoofVeg';'WBRoofTot'};

for i=1:size(WBRoofNames,1)
	WBRoof.(cell2mat(WBRoofNames(i)))	=	zeros(n,1,m);
end


WBCanyonIndvNames	=	{'WB_In_tree';'WB_In_gveg';'WB_In_gimp';'WB_In_gbare';...
						'WB_Pond_gveg';'WB_Soil_gimp';'WB_Soil_gbare';'WB_Soil_gveg'};

for i=1:size(WBCanyonIndvNames,1)
	WBCanyonIndv.(cell2mat(WBCanyonIndvNames(i)))	=	zeros(n,1,m);
end


WBCanyonTotNames	=	{'WBsurf_tree';'WBsurf_imp';'WBsurf_bare';'WBsurf_veg';...
						'WBsoil_imp';'WBsoil_bare';'WBsoil_veg';'WBimp_tot';'WBbare_tot';'WBveg_tot';...
						'WBcanyon_flux';'WBtree_level';'WBground_level';'WBsoil_level';'WBcanyon_level'};

for i=1:size(WBCanyonTotNames,1)
	WBCanyonTot.(cell2mat(WBCanyonTotNames(i)))	=	zeros(n,1,m);
end


% Energy Balance
% EBGroundImp	=	EBalance_groundimp
% EBGroundBare	=	EBalance_groundbare
% EBGroundVeg	=	EBalance_groundveg
% EBTree		=	EBalance_tree
% EBWallSun		=	EBalance_wallsun
% EBWallShade	=	EBalance_wallshade
% EBWallSunInt	=	EBalance_wallsun_interior
% EBWallShadeInt=	EBalance_wallshade_interior
% EBCanyonT		=	EBalance_canyon_temp
% EBCanyonQ		=	EBalance_canyon_humid

EBNames	=	{'EBRoofImp';'EBRoofVeg';'EBGroundImp';'EBGroundBare';...
					'EBGroundVeg';'EBTree';'EBWallSun';'EBWallShade';'EBWallSunInt';...
					'EBWallShadeInt';'EBCanyonT';'EBCanyonQ'};

for i=1:size(EBNames,1)
	EB.(cell2mat(EBNames(i)))	=	zeros(n,1,m);
end

%% Wind speed
WindNames	=	{'u_Hcan';'u_Zref_und'};

for i=1:size(WindNames,1)
	Wind.(cell2mat(WindNames(i)))	=	zeros(n,1,m);
end

%% Success of energy balance solver
Solver				=	[];
Solver.SuccessRoof	=	zeros(n,1,m);
Solver.SuccessCanyon=	zeros(n,1,m);
Solver.ValuesRoof	=	zeros(n,4,m);
Solver.ValuesCanyon	=	zeros(n,10,m);
Solver.TsolverRoof	=	zeros(n,4,m);
Solver.TsolverCanyon=	zeros(n,10,m);

%% Temperature and humidity at 2m canyon height
Results2m		=	[];
Results2m.T2m	=	zeros(n,1,m);
Results2m.q2m	=	zeros(n,1,m);
Results2m.e_T2m	=	zeros(n,1,m);
Results2m.RH_T2m=	zeros(n,1,m);
Results2m.qcan	=	zeros(n,1,m);
Results2m.e_Tcan=	zeros(n,1,m);
Results2m.RH_Tcan=	zeros(n,1,m);


%% Start calculation
tic
Opt_Solv = optimoptions('lsqnonlin','Display','off');

for ittm = 1:m

%% Initialize variables
[Gemeotry_m,ParTree,geometry,FractionsRoof,FractionsGround,...
	WallLayers,ParSoilRoof,ParSoilGround,ParInterceptionTree,...
	PropOpticalRoof,PropOpticalGround,PropOpticalWall,PropOpticalTree,...
	ParThermalRoof,ParThermalGround,ParThermalWall,ParThermalTree,...
	ParVegRoof,ParVegGround,ParVegTree]=feval(strcat('data_functions.Data_UEHM_site_',Name_Site),ittm);

% [Gemeotry_m,ParTree,geometry,FractionsRoof,FractionsGround,...
% 	WallLayers,ParSoilRoof,ParSoilGround,ParInterceptionTree,...
% 	PropOpticalRoof,PropOpticalGround,PropOpticalWall,PropOpticalTree,...
% 	ParThermalRoof,ParThermalGround,ParThermalWall,ParThermalTree,...
% 	ParVegRoof,ParVegGround,ParVegTree]=data_functions.Data_UEHM_site(ittm);

%% Calculate view factors
[ViewFactor]=ray_tracing.VFUrbanCanyon(OPTION_RAY,Name_Site,Gemeotry_m,geometry);

if FractionsRoof.fimp==1
Vwater.VRoofSoilVeg(1,:,ittm)	=	repmat(0.*ParSoil.Roof.dz,1,1,1);		% Starting with dry soil
Owater.OwRoofSoilVeg(1,:,ittm)	=	0;				% Starting with dry soil
end

if FractionsGround.fimp==1
Vwater.VGroundSoilImp(1,:,ittm)		=	repmat(0.*ParSoil.Ground.dz,1,1,1);	% Starting with dry soil
Owater.OwGroundSoilImp(1,:,ittm)	=	0;				% Starting with dry soil
Vwater.VGroundSoilBare(1,:,ittm)	=	repmat(0.*ParSoil.Ground.dz,1,1,1);	% Starting with dry soil
Owater.OwGroundSoilBare(1,:,ittm)	=	0;				% Starting with dry soil
Vwater.VGroundSoilVeg(1,:,ittm)		=	repmat(0.*ParSoil.Ground.dz,1,1,1);	% Starting with dry soil
Owater.OwGroundSoilVeg(1,:,ittm)	=	0;				% Starting with dry soil
end

for ittn	= 1:n
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[SunPosition,MeteoData,HumidityAtm,Anthropogenic,location,ParCalculation]...
	=feval(strcat('data_functions.UEHMForcingData_',Name_Site),MeteoDataRaw,ittn);

% [SunPosition,MeteoData,HumidityAtm,Anthropogenic,location,ParCalculation]...
% 	=data_functions.UEHMForcingData(MeteoDataRaw,ittn);

	for i=1:size(TempVecNames,1)
		if ittn==1
			TempVec_ittm.(cell2mat(TempVecNames(i)))	=	TempVec.(cell2mat(TempVecNames(i)))(1,:,ittm); 
		else
			TempVec_ittm.(cell2mat(TempVecNames(i)))	=	TempVec.(cell2mat(TempVecNames(i)))(ittn-1,:,ittm);
		end
	end
	
	for i=1:size(IntNames,1)
		if ittn==1
			Int_ittm.(cell2mat(IntNames(i)))	=	Int.(cell2mat(IntNames(i)))(1,:,ittm); 
		else
			Int_ittm.(cell2mat(IntNames(i)))	=	Int.(cell2mat(IntNames(i)))(ittn-1,:,ittm); 
		end
	end
	
	for i=1:size(ExWaterNames,1)
		if ittn==1
			ExWater_ittm.(cell2mat(ExWaterNames(i)))	=	ExWater.(cell2mat(ExWaterNames(i)))(1,:,ittm);
		else
			ExWater_ittm.(cell2mat(ExWaterNames(i)))	=	ExWater.(cell2mat(ExWaterNames(i)))(ittn-1,:,ittm);
		end
	end
	
	VwaterNames	=	fieldnames(Vwater);
	for i=1:size(VwaterNames,1)
		if ittn==1
			Vwater_ittm.(cell2mat(VwaterNames(i)))	=	Vwater.(cell2mat(VwaterNames(i)))(1,:,ittm);
		else
			Vwater_ittm.(cell2mat(VwaterNames(i)))	=	Vwater.(cell2mat(VwaterNames(i)))(ittn-1,:,ittm);
		end
	end
	
	OwaterNames	=	fieldnames(Owater);
	for i=1:size(OwaterNames,1)
		if ittn==1
			Owater_ittm.(cell2mat(OwaterNames(i)))	=	Owater.(cell2mat(OwaterNames(i)))(1,:,ittm);
		else
			Owater_ittm.(cell2mat(OwaterNames(i)))	=	Owater.(cell2mat(OwaterNames(i)))(ittn-1,:,ittm);
		end
	end
	
	for i=1:size(SoilPotWNames,1)
		if ittn==1
			SoilPotW_ittm.(cell2mat(SoilPotWNames(i)))	=	SoilPotW.(cell2mat(SoilPotWNames(i)))(1,:,ittm);
		else
			SoilPotW_ittm.(cell2mat(SoilPotWNames(i)))	=	SoilPotW.(cell2mat(SoilPotWNames(i)))(ittn-1,:,ittm);
		end
	end
	
	for i=1:size(CiCO2LeafNames,1)
		if ittn==1
			CiCO2Leaf_ittm.(cell2mat(CiCO2LeafNames(i)))	=	CiCO2Leaf.(cell2mat(CiCO2LeafNames(i)))(1,:,ittm); 
		else
			CiCO2Leaf_ittm.(cell2mat(CiCO2LeafNames(i)))	=	CiCO2Leaf.(cell2mat(CiCO2LeafNames(i)))(ittn-1,:,ittm); 
		end
	end
	
	for i=1:size(HumidityNames,1)
		if ittn==1
			Humidity_ittm.(cell2mat(HumidityNames(i)))	=	Humidity.(cell2mat(HumidityNames(i)))(1,:,ittm); 
		else
			Humidity_ittm.(cell2mat(HumidityNames(i)))	=	Humidity.(cell2mat(HumidityNames(i)))(ittn-1,:,ittm); 
		end
	end
	
	for i=1:size(TempDampNames,1)
		if ittn==1
			TempDamp_ittm.(cell2mat(TempDampNames(i)))	=	TempDamp.(cell2mat(TempDampNames(i)))(1,:,ittm);
		else
			TempDamp_ittm.(cell2mat(TempDampNames(i)))	=	TempDamp.(cell2mat(TempDampNames(i)))(ittn-1,:,ittm);
		end
	end
	
	for i=1:size(RunonNames,1)
		if ittn==1
			Runon_ittm.(cell2mat(RunonNames(i)))	=	Runon.(cell2mat(RunonNames(i)))(1,:,ittm); 
		else
			Runon_ittm.(cell2mat(RunonNames(i)))	=	Runon.(cell2mat(RunonNames(i)))(ittn-1,:,ittm); 
		end
	end
	
	QinlatNames	=	fieldnames(Qinlat);
	for i=1:size(QinlatNames,1)
		if ittn==1
			Qinlat_ittm.(cell2mat(QinlatNames(i)))	=	Qinlat.(cell2mat(QinlatNames(i)))(1,:,ittm); 
		else
			Qinlat_ittm.(cell2mat(QinlatNames(i)))	=	Qinlat.(cell2mat(QinlatNames(i)))(ittn-1,:,ittm); 
		end
	end
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
% Calculate Energy Budget of Roof
	ittn_only	=	1;
	[TR,fvalR,exitflagR]=fSolver_roof(TempVec_ittm,MeteoData,ittn_only,...
	Int_ittm,ExWater_ittm,Vwater_ittm,Owater_ittm,SoilPotW_ittm,CiCO2Leaf_ittm,Opt_Solv,...
	Gemeotry_m,FractionsRoof,ParSoilRoof,PropOpticalRoof,ParThermalRoof,ParVegRoof,...
	HumidityAtm,Anthropogenic,ParCalculation);

	TempVec.TRoofImp(ittn,1,ittm)		=	TR(1,1);
	TempVec.TRoofVeg(ittn,1,ittm)		=	TR(1,2);
	TempVec.TRoofIntImp(ittn,1,ittm)	=	TR(1,3);
	TempVec.TRoofIntVeg(ittn,1,ittm)	=	TR(1,4);
	
% Calculate Energy Budget of the Canyon
	[TC,fvalC,exitflagC]=fSolver_canyon(TempVec_ittm,Humidity_ittm,MeteoData,ittn_only,...
	Int_ittm,ExWater_ittm,Vwater_ittm,Owater_ittm,SoilPotW_ittm,CiCO2Leaf_ittm,TempDamp_ittm,ViewFactor,Opt_Solv,...
	Gemeotry_m,ParTree,geometry,FractionsGround,...
	WallLayers,ParSoilGround,ParInterceptionTree,...
	PropOpticalGround,PropOpticalWall,PropOpticalTree,...
	ParThermalGround,ParThermalWall,ParVegGround,ParVegTree,...
	SunPosition,HumidityAtm,Anthropogenic,ParCalculation);

	TempVec.TGroundImp(ittn,1,ittm)		=	TC(1,1);
	TempVec.TGroundBare(ittn,1,ittm)	=	TC(1,2);
	TempVec.TGroundVeg(ittn,1,ittm)		=	TC(1,3);
	TempVec.TTree(ittn,1,ittm)			=	TC(1,6);
	TempVec.TWallSun(ittn,1,ittm)		=	TC(1,4);
	TempVec.TWallShade(ittn,1,ittm)		=	TC(1,5);
	TempVec.TWallIntSun(ittn,1,ittm)	=	TC(1,7);
	TempVec.TWallIntShade(ittn,1,ittm)	=	TC(1,8);
	TempVec.TCanyon(ittn,1,ittm)		=	TC(1,9);
	Humidity.CanyonSpecific(ittn,1,ittm)=	TC(1,10);

% Calculate Energy and Water fluxes Roof
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[SWRabsRoofImp,SWRabsRoofVeg,SWRabsTotalRoof,...
		SWRoutRoofImp,SWRoutRoofVeg,SWRoutTotalRoof,...
		SWRinRoofImp,SWRinRoofVeg,SWRinTotalRoof,...
		SWREBRoofImp,SWREBRoofVeg,SWREBTotalRoof,...
		LWRabsRoofVeg,LWRabsRoofImp,LWRabsTotalRoof,...
		LWRoutRoofVeg,LWRoutRoofImp,LWRoutTotalRoof,...
		LWRinRoofImp,LWRinRoofVeg,LWRinTotalRoof,...
		LWREBRoofImp,LWREBRoofVeg,LWREBTotalRoof,...
		HfluxRoofImp,HfluxRoofVeg,HfluxRoof,...
		LEfluxRoofImp,LEfluxRoofVegInt,LEfluxRoofVegPond,...
		LEfluxRoofVegSoil,LTEfluxRoofVeg,LEfluxRoofVeg,LEfluxRoof,...
		G1RoofImp,G2RoofImp,dsRoofImp,G1RoofVeg,G2RoofVeg,dsRoofVeg,G1Roof,G2Roof,dsRoof,...
		raRooftoAtm,rb_LRoof,rap_LRoof,r_soilRoof,rs_sunRoof,rs_shdRoof,...
		EfluxRoofImp,EfluxRoofVegInt,EfluxRoofVegPond,...
		EfluxRoofVegSoil,TEfluxRoofVeg,EfluxRoofVeg,EfluxRoof,...
		... % Water fluxes
		QRoofImp,QRoofVegDrip,QRoofVegPond,LkRoofImp,LkRoofVeg,LkRoof,QRoofVegSoil,RunoffRoofTot,RunonRoofTot,...
		IntRoofImp,IntRoofVegPlant,IntRoofVegGround,dInt_dtRoofImp,dInt_dtRoofVegPlant,dInt_dtRoofVegGround,...
		IntRooftot,dInt_dtRooftot,dVRoofSoilVeg_dt,...
		fRoofVeg,VRoofSoilVeg,OwRoofSoilVeg,OSwRoofSoilVeg,ExWaterRoofVeg_H,SoilPotWRoofVeg_H,SoilPotWRoofVeg_L,ExWaterRoofVeg_L,...
		CiCO2LeafRoofVegSun,CiCO2LeafRoofVegShd,...
		WBRoofVegInVeg,WBRoofVegInGround,WBRoofVegSoil,...
		...
		EBRoofImp,EBRoofVeg,Yroof,WBRoofImp,WBRoofVeg,WBRoofTot]...
		=EB_WB_roof(TR,TempVec_ittm,MeteoData,ittn_only,...
		Int_ittm,ExWater_ittm,Vwater_ittm,Owater_ittm,SoilPotW_ittm,CiCO2Leaf_ittm,Runon_ittm,...
		Gemeotry_m,FractionsRoof,ParSoilRoof,PropOpticalRoof,ParThermalRoof,ParVegRoof,...
		HumidityAtm,Anthropogenic,ParCalculation);
	
% Calculate Energy and Water fluxes Canyon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[SWRin_t,SWRout_t,SWRabs_t,SWRabsDir_t,SWRabsDiff_t,SWREB_t,...
		LWRin_t,LWRout_t,LWRabs_t,LWREB_t,...
		HfluxGroundImp,HfluxGroundBare,HfluxGroundVeg,HfluxTree,HfluxGround,...
		EfluxGroundImp,EfluxGroundBarePond,EfluxGroundBareSoil,EfluxGroundVegInt,...
		EfluxGroundVegPond,EfluxGroundVegSoil,TEfluxGroundVeg,EfluxTreeInt,TEfluxTree,...
		EfluxGroundBare,EfluxGroundVeg,EfluxGround,EfluxTree,...
		LEfluxGroundImp,LEfluxGroundBarePond,LEfluxGroundBareSoil,LEfluxGroundVegInt,...
		LEfluxGroundVegPond,LEfluxGroundVegSoil,LTEfluxGroundVeg,LEfluxTreeInt,LTEfluxTree,...
		LEfluxGroundBare,LEfluxGroundVeg,LEfluxGround,LEfluxTree,...
		CiCO2LeafTreeSun,CiCO2LeafTreeShd,CiCO2LeafGroundVegSun,CiCO2LeafGroundVegShd,...
		raCanyontoAtm,rap_can,rap_Htree_In,rb_HGround,rb_LGround,...
		r_soilGroundbare,r_soilGroundveg,alp_soilGroundbare,alp_soilGroundveg,...
		rs_sunGround,rs_shdGround,rs_sunTree,rs_shdTree,...
		Fsun_L,Fshd_L,dw_L,RES_w1,RES_w2,rap_W1_In,rap_W2_In,...
		HfluxWallSun,HfluxWallShade,EfluxWallSun,EfluxWallShade,LEfluxWallSun,LEfluxWallShade,HfluxCanyon,LEfluxCanyon,EfluxCanyon,...
		G1WallSun,G2WallSun,dsWallSun,G1WallShade,G2WallShade,dsWallShade,...
		G1GroundImp,TDampGroundImp,G1GroundBare,TDampGroundBare,G1GroundVeg,TDampGroundVeg,GTree,TDampTree,G1Ground,G1Canyon,G2Canyon,...
		dsGroundImp,dsGroundBare,dsGroundVeg,dsTree,dsCanyonAir,Ycanyon,...
		... %Water fluxes
		QTree,IntTree,dInt_dtTree,QGroundVegDrip,IntGroundVegPlant,dInt_dtGroundVegPlant,...
		QGroundImp,IntGroundImp,dInt_dtGroundImp,fGroundImp,QGroundBarePond,IntGroundBare,dInt_dtGroundBare,fGroundBare,...
		QGroundVegPond,IntGroundVegGround,dInt_dtGroundVegGround,fGroundVeg,...
		...
		VGroundSoilImp,OwGroundSoilImp,OSwGroundSoilImp,LkGroundImp,SoilPotWGroundImp_H,SoilPotWGroundImp_L,...
		ExWaterGroundImp_H,ExWaterGroundImp_L,Rd_gimp,TEgveg_imp,TEtree_imp,...
		Egimp_soil,dVGroundSoilImp_dt,Psi_Soil_gimp,Kf_gimp,...
		...
		VGroundSoilBare,OwGroundSoilBare,OSwGroundSoilBare,LkGroundBare,SoilPotWGroundBare_H,SoilPotWGroundBare_L,...
		ExWaterGroundBare_H,ExWaterGroundBare_L,QGroundBareSoil,TEgveg_bare,TEtree_bare,...
		Egbare_Soil,dVGroundSoilBare_dt,Psi_soil_gbare,Kf_gbare,...
		...
		VGroundSoilVeg,OwGroundSoilVeg,OSwGroundSoilVeg,LkGroundVeg,SoilPotWGroundVeg_H,SoilPotWGroundVeg_L,...
		ExWaterGroundVeg_H,ExWaterGroundVeg_L,QGroundVegSoil,TEgveg_veg,TEtree_veg,...
		Egveg_Soil,dVGroundSoilVeg_dt,Psi_soil_gveg,Kf_gveg,...
		...
		Qin_imp,Qin_bare,Qin_veg,Qin_bare2imp,Qin_bare2veg,Qin_imp2bare,Qin_imp2veg,Qin_veg2imp,Qin_veg2bare,...
		...
		VGroundSoilTot,OwGroundSoilTot,OSwGroundSoilTot,LkGround,Rd,dVGroundSoilTot_dt,SoilPotWGroundTot_L,ExWaterGroundTot_L,TEgveg_tot,SoilPotWGroundTot_H,ExWaterGroundTot_H,...
		TEtree_tot,EB_TEtree,EB_TEgveg,WBIndv,WBTot,...
		RunoffGroundTot,RunonGroundTot,Etot,DeepGLk,StorageTot,...
		...
		EBGroundImp,EBGroundBare,EBGroundVeg,EBTree,EBWallSun,EBWallShade,EBWallSunInt,EBWallShadeInt,EBCanyonT,EBCanyonQ,...
		HumidityCan,HumidityAtm,u_Hcan,u_Zref_und,T2m,q2m,e_T2m,RH_T2m,qcan,e_Tcan,RH_Tcan]...
		=EB_WB_canyon(TC,TempVec_ittm,MeteoData,ittn_only,...
		Int_ittm,ExWater_ittm,Vwater_ittm,Owater_ittm,SoilPotW_ittm,...
		CiCO2Leaf_ittm,TempDamp_ittm,Runon_ittm,Qinlat_ittm,ViewFactor,...
		Gemeotry_m,ParTree,geometry,FractionsGround,...
		WallLayers,ParSoilGround,ParInterceptionTree,...
		PropOpticalGround,PropOpticalWall,PropOpticalTree,...
		ParThermalGround,ParThermalWall,ParVegGround,ParVegTree,...
		SunPosition,HumidityAtm,Anthropogenic,ParCalculation);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SWRabs_t.SWRabsTotalUrban	=	geometry.wroof_norm*SWRabsTotalRoof + geometry.wcanyon_norm*SWRabs_t.SWRabsTotalCanyon;
SWRin_t.SWRinTotalUrban		=	geometry.wroof_norm*SWRinTotalRoof + geometry.wcanyon_norm*SWRin_t.SWRinTotalCanyon;
SWRout_t.SWRoutTotalUrban	=	geometry.wroof_norm*SWRoutTotalRoof + geometry.wcanyon_norm*SWRout_t.SWRoutTotalCanyon;
SWREB_t.SWREBTotalUrban		=	geometry.wroof_norm*SWREBTotalRoof + geometry.wcanyon_norm*SWREB_t.SWREBTotalCanyon;

LWRabs_t.LWRabsTotalUrban	=	geometry.wroof_norm*LWRabsTotalRoof + geometry.wcanyon_norm*LWRabs_t.LWRabsTotalCanyon;
LWRin_t.LWRinTotalUrban		=	geometry.wroof_norm*LWRinTotalRoof + geometry.wcanyon_norm*LWRin_t.LWRinTotalCanyon;
LWRout_t.LWRoutTotalUrban	=	geometry.wroof_norm*LWRoutTotalRoof + geometry.wcanyon_norm*LWRout_t.LWRoutTotalCanyon;
LWREB_t.LWREBTotalUrban		=	geometry.wroof_norm*LWREBTotalRoof + geometry.wcanyon_norm*LWREB_t.LWREBTotalCanyon;

HfluxUrban	=	geometry.wroof_norm*HfluxRoof + geometry.wcanyon_norm*HfluxCanyon;
LEfluxUrban	=	geometry.wroof_norm*LEfluxRoof + geometry.wcanyon_norm*LEfluxCanyon;
G1Urban		=	geometry.wroof_norm*G1Roof + geometry.wcanyon_norm*G1Canyon;
G2Urban		=	geometry.wroof_norm*G2Roof + geometry.wcanyon_norm*G2Canyon;
EfluxUrban	=	geometry.wroof_norm*EfluxRoof + geometry.wcanyon_norm*EfluxCanyon;
RunonUrban	=	geometry.wroof_norm*RunonRoofTot + geometry.wcanyon_norm*RunonGroundTot;
RunoffUrban	=	geometry.wroof_norm*RunoffRoofTot + geometry.wcanyon_norm*RunoffGroundTot;
LkUrban		=	geometry.wroof_norm*LkRoof + geometry.wcanyon_norm*LkGround;
	
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tracking the success of the energy balance solver
Solver.SuccessRoof(ittn,:,ittm)		=	exitflagR;
Solver.SuccessCanyon(ittn,:,ittm)	=	exitflagC;
Solver.ValuesRoof(ittn,:,ittm)		=	fvalR;
Solver.ValuesCanyon(ittn,:,ittm)	=	fvalC;
Solver.TsolverRoof(ittn,:,ittm)		=	TR;
Solver.TsolverCanyon(ittn,:,ittm)	=	TC;

% Results 2m
Results2m.T2m(ittn,:,ittm)		=	T2m;
Results2m.q2m(ittn,:,ittm)		=	q2m;
Results2m.e_T2m(ittn,:,ittm)	=	e_T2m;
Results2m.RH_T2m(ittn,:,ittm)	=	RH_T2m;
Results2m.qcan(ittn,:,ittm)		=	qcan;
Results2m.e_Tcan(ittn,:,ittm)	=	e_Tcan;
Results2m.RH_Tcan(ittn,:,ittm)	=	RH_Tcan;

% Anthropogenic time series
Anthropo.Tb(ittn,:,ittm)				=	Anthropogenic.Tb;
Anthropo.Qf_canyon(ittn,:,ittm)			=	Anthropogenic.Qf_canyon;
Anthropo.Qf_roof(ittn,:,ittm)			=	Anthropogenic.Qf_roof;
Anthropo.Waterf_canyonVeg(ittn,:,ittm)	=	Anthropogenic.Waterf_canyonVeg;
Anthropo.Waterf_canyonBare(ittn,:,ittm)	=	Anthropogenic.Waterf_canyonBare;
Anthropo.Waterf_roof(ittn,:,ittm)		=	Anthropogenic.Waterf_roof;


% Humidity
for i=1:6
	Humidity.(cell2mat(HumidityNames(i)))(ittn,1,ittm)=	HumidityCan.(cell2mat(HumidityNames(i)));
end

for i=7:12
	Humidity.(cell2mat(HumidityNames(i)))(ittn,1,ittm)=	HumidityAtm.(cell2mat(HumidityNames(i)));
end

% Shortwave radiation
for i=1:3
	SWRabs.(cell2mat(SWRabsNames(i)))(ittn,1,ittm)	=	eval(cell2mat(SWRabsNames(i)));
	SWRin.(cell2mat(SWRinNames(i)))(ittn,1,ittm)	=	eval(cell2mat(SWRinNames(i)));
	SWRout.(cell2mat(SWRoutNames(i)))(ittn,1,ittm)	=	eval(cell2mat(SWRoutNames(i)));
	SWREB.(cell2mat(SWREBNames(i)))(ittn,1,ittm)	=	eval(cell2mat(SWREBNames(i)));
end
for i=4:size(SWRabsNames,1)
	SWRabs.(cell2mat(SWRabsNames(i)))(ittn,1,ittm)	=	SWRabs_t.(cell2mat(SWRabsNames(i)))(1,1);
	SWRin.(cell2mat(SWRinNames(i)))(ittn,1,ittm)	=	SWRin_t.(cell2mat(SWRinNames(i)))(1,1);
	SWRout.(cell2mat(SWRoutNames(i)))(ittn,1,ittm)	=	SWRout_t.(cell2mat(SWRoutNames(i)))(1,1);
	SWREB.(cell2mat(SWREBNames(i)))(ittn,1,ittm)	=	SWREB_t.(cell2mat(SWREBNames(i)))(1,1);
end

% Longwave radiation
for i=1:3
	LWRabs.(cell2mat(LWRabsNames(i)))(ittn,1,ittm)	=	eval(cell2mat(LWRabsNames(i)));
	LWRin.(cell2mat(LWRinNames(i)))(ittn,1,ittm)	=	eval(cell2mat(LWRinNames(i)));
	LWRout.(cell2mat(LWRoutNames(i)))(ittn,1,ittm)	=	eval(cell2mat(LWRoutNames(i)));
	LWREB.(cell2mat(LWREBNames(i)))(ittn,1,ittm)	=	eval(cell2mat(LWREBNames(i)));
end
for i=4:size(LWRabsNames,1)
	LWRabs.(cell2mat(LWRabsNames(i)))(ittn,1,ittm)	=	LWRabs_t.(cell2mat(LWRabsNames(i)))(1,1);
	LWRin.(cell2mat(LWRinNames(i)))(ittn,1,ittm)	=	LWRin_t.(cell2mat(LWRinNames(i)))(1,1);
	LWRout.(cell2mat(LWRoutNames(i)))(ittn,1,ittm)	=	LWRout_t.(cell2mat(LWRoutNames(i)))(1,1);
	LWREB.(cell2mat(LWREBNames(i)))(ittn,1,ittm)	=	LWREB_t.(cell2mat(LWREBNames(i)))(1,1);
end

% Sensible heat
for i=1:size(HfluxNames,1)
	Hflux.(cell2mat(HfluxNames(i)))(ittn,1,ittm)	=	eval(cell2mat(HfluxNames(i)));
end

% Latent heat
for i=1:size(LEfluxNames,1)
	LEflux.(cell2mat(LEfluxNames(i)))(ittn,1,ittm)	=	eval(cell2mat(LEfluxNames(i)));
end

% Ground heat flux
for i=1:size(GfluxNames,1)
	Gflux.(cell2mat(GfluxNames(i)))(ittn,1,ittm)	=	eval(cell2mat(GfluxNames(i)));
end

% Heat Storage
for i=1:size(dStorageNames,1)
	dStorage.(cell2mat(dStorageNames(i)))(ittn,1,ittm)	=	eval(cell2mat(dStorageNames(i)));
end

% Dampening temperature
for i=1:size(TempDampNames,1)
	TempDamp.(cell2mat(TempDampNames(i)))(ittn,1,ittm)	=	eval(cell2mat(TempDampNames(i)));
end

for i=1:size(RESNames,1)
	RES.(cell2mat(RESNames(i)))(ittn,1,ittm)	=	eval(cell2mat(RESNames(i)));
end

for i=1:size(EfluxNames,1)
	Eflux.(cell2mat(EfluxNames(i)))(ittn,1,ittm)	=	eval(cell2mat(EfluxNames(i)));
end

for i=1:size(RunoffNames,1)
	Runoff.(cell2mat(RunoffNames(i)))(ittn,1,ittm)	=	eval(cell2mat(RunoffNames(i)));
end

for i=1:size(RunonNames,1)
	Runon.(cell2mat(RunonNames(i)))(ittn,1,ittm)	=	eval(cell2mat(RunonNames(i)));
end

for i=1:size(LeakageNames,1)
	Leakage.(cell2mat(LeakageNames(i)))(ittn,1,ittm)	=	eval(cell2mat(LeakageNames(i)));
end

for i=1:size(IntNames,1)
	Int.(cell2mat(IntNames(i)))(ittn,1,ittm)	=	eval(cell2mat(IntNames(i)));
end

for i=1:size(dInt_dtNames,1)
	dInt_dt.(cell2mat(dInt_dtNames(i)))(ittn,1,ittm)	=	eval(cell2mat(dInt_dtNames(i)));
end

for i=1:size(InfiltrationNames,1)
	Infiltration.(cell2mat(InfiltrationNames(i)))(ittn,1,ittm)=	eval(cell2mat(InfiltrationNames(i)));
end

Vwater.VRoofSoilVeg(ittn,:,ittm)	=	VRoofSoilVeg;
Vwater.VGroundSoilImp(ittn,:,ittm)	=	VGroundSoilImp;
Vwater.VGroundSoilBare(ittn,:,ittm)	=	VGroundSoilBare;
Vwater.VGroundSoilVeg(ittn,:,ittm)	=	VGroundSoilVeg;
Vwater.VGroundSoilTot(ittn,:,ittm)	=	VGroundSoilTot;

dVwater_dt.dVRoofSoilVeg_dt(ittn,:,ittm)	=	dVRoofSoilVeg_dt;
dVwater_dt.dVGroundSoilImp_dt(ittn,:,ittm)	=	dVGroundSoilImp_dt;
dVwater_dt.dVGroundSoilBare_dt(ittn,:,ittm)	=	dVGroundSoilBare_dt;
dVwater_dt.dVGroundSoilVeg_dt(ittn,:,ittm)	=	dVGroundSoilVeg_dt;
dVwater_dt.dVGroundSoilTot_dt(ittn,:,ittm)	=	dVGroundSoilTot_dt;

Owater.OwRoofSoilVeg(ittn,:,ittm)	=	OwRoofSoilVeg;
Owater.OwGroundSoilImp(ittn,:,ittm)	=	OwGroundSoilImp;
Owater.OwGroundSoilBare(ittn,:,ittm)=	OwGroundSoilBare;
Owater.OwGroundSoilVeg(ittn,:,ittm)	=	OwGroundSoilVeg;
Owater.OwGroundSoilTot(ittn,:,ittm)	=	OwGroundSoilTot;

OSwater.OSwRoofSoilVeg(ittn,:,ittm)	=	OSwRoofSoilVeg;
OSwater.OSwGroundSoilImp(ittn,:,ittm)=	OSwGroundSoilImp;
OSwater.OSwGroundSoilBare(ittn,:,ittm)=	OSwGroundSoilBare;
OSwater.OSwGroundSoilVeg(ittn,:,ittm)=	OSwGroundSoilVeg;
OSwater.OSwGroundSoilTot(ittn,:,ittm)=	OSwGroundSoilTot;

% Lateral soil water flux
Qinlat.Qin_bare2imp(ittn,:,ittm)=	Qin_bare2imp;
Qinlat.Qin_veg2imp(ittn,:,ittm)	=	Qin_veg2imp;
Qinlat.Qin_veg2bare(ittn,:,ittm)=	Qin_veg2bare;
Qinlat.Qin_imp2bare(ittn,:,ittm)=	Qin_imp2bare;
Qinlat.Qin_bare2veg(ittn,:,ittm)=	Qin_bare2veg;
Qinlat.Qin_imp2veg(ittn,:,ittm)	=	Qin_imp2veg;
Qinlat.Qin_imp(ittn,:,ittm)		=	Qin_imp;
Qinlat.Qin_bare(ittn,:,ittm)	=	Qin_bare;
Qinlat.Qin_veg(ittn,:,ittm)		=	Qin_veg;


for i=1:size(ExWaterNames,1)
	ExWater.(cell2mat(ExWaterNames(i)))(ittn,:,ittm)=	eval(cell2mat(ExWaterNames(i)));
end

for i=1:size(SoilPotWNames,1)
	SoilPotW.(cell2mat(SoilPotWNames(i)))(ittn,1,ittm)=	eval(cell2mat(SoilPotWNames(i)));
end

for i=1:size(CiCO2LeafNames,1)
	CiCO2Leaf.(cell2mat(CiCO2LeafNames(i)))(ittn,1,ittm)=	eval(cell2mat(CiCO2LeafNames(i)));
end

for i=1:size(EBNames,1)
	EB.(cell2mat(EBNames(i)))(ittn,1,ittm)=	eval(cell2mat(EBNames(i)));
end

WBRoofVegSoil	=	sum(WBRoofVegSoil);
for i=1:size(WBRoofNames,1)
	WBRoof.(cell2mat(WBRoofNames(i)))(ittn,:,ittm)=	eval(cell2mat(WBRoofNames(i)));
end

WBCanyonIndv.WB_In_tree(ittn,:,ittm)	=	WBIndv.WB_In_tree;
WBCanyonIndv.WB_In_gveg(ittn,:,ittm)	=	WBIndv.WB_In_gveg;
WBCanyonIndv.WB_In_gimp(ittn,:,ittm)	=	WBIndv.WB_In_gimp;
WBCanyonIndv.WB_In_gbare(ittn,:,ittm)	=	WBIndv.WB_In_gbare;
WBCanyonIndv.WB_Pond_gveg(ittn,:,ittm)	=	WBIndv.WB_Pond_gveg;
WBCanyonIndv.WB_Soil_gimp(ittn,:,ittm)	=	nansum(WBIndv.WB_Soil_gimp);
WBCanyonIndv.WB_Soil_gbare(ittn,:,ittm)	=	sum(WBIndv.WB_Soil_gbare);
WBCanyonIndv.WB_Soil_gveg(ittn,:,ittm)	=	sum(WBIndv.WB_Soil_gveg);

WBCanyonTot.WBsurf_tree(ittn,:,ittm)	=	WBTot.WBsurf_tree;
WBCanyonTot.WBsurf_imp(ittn,:,ittm)		=	WBTot.WBsurf_imp;
WBCanyonTot.WBsurf_bare(ittn,:,ittm)	=	WBTot.WBsurf_bare;
WBCanyonTot.WBsurf_veg(ittn,:,ittm)		=	WBTot.WBsurf_veg;
WBCanyonTot.WBsoil_imp(ittn,:,ittm)		=	WBTot.WBsoil_imp;
WBCanyonTot.WBsoil_bare(ittn,:,ittm)	=	WBTot.WBsoil_bare;
WBCanyonTot.WBsoil_veg(ittn,:,ittm)		=	WBTot.WBsoil_veg;
WBCanyonTot.WBimp_tot(ittn,:,ittm)		=	WBTot.WBimp_tot;
WBCanyonTot.WBbare_tot(ittn,:,ittm)		=	WBTot.WBbare_tot;
WBCanyonTot.WBveg_tot(ittn,:,ittm)		=	WBTot.WBveg_tot;
WBCanyonTot.WBcanyon_flux(ittn,:,ittm)	=	WBTot.WBcanyon_flux;
WBCanyonTot.WBtree_level(ittn,:,ittm)	=	WBTot.WBtree_level;
WBCanyonTot.WBground_level(ittn,:,ittm)	=	WBTot.WBground_level;
WBCanyonTot.WBsoil_level(ittn,:,ittm)	=	WBTot.WBsoil_level;
WBCanyonTot.WBcanyon_level(ittn,:,ittm)	=	WBTot.WBcanyon_level;

for i=1:size(EBNames,1)
	EB.(cell2mat(EBNames(i)))(ittn,1,ittm)=	eval(cell2mat(EBNames(i)));
end

for i=1:size(WindNames,1)
	Wind.(cell2mat(WindNames(i)))(ittn,1,ittm)=	eval(cell2mat(WindNames(i)));
end


if mod(ittn,10)==0
disp(strcat('iter=',num2str(ittn),' iterm=',num2str(ittm)));
end 

end

Gemeotry_m_Out(ittm)			=	Gemeotry_m;
ParTree_Out(ittm)				=	ParTree;
geometry_Out(ittm)				=	geometry;
FractionsRoof_Out(ittm)			=	FractionsRoof;
FractionsGround_Out(ittm)		=	FractionsGround;
WallLayers_Out(ittm)			=	WallLayers;
ParSoilRoof_Out(ittm)			=	ParSoilRoof;
ParSoilGround_Out(ittm)			=	ParSoilGround;
ParInterceptionTree_Out(ittm)	=	ParInterceptionTree;
PropOpticalRoof_Out(ittm)		=	PropOpticalRoof;
PropOpticalGround_Out(ittm)		=	PropOpticalGround;
PropOpticalWall_Out(ittm)		=	PropOpticalWall;
PropOpticalTree_Out(ittm)		=	PropOpticalTree;
ParThermalRoof_Out(ittm)		=	ParThermalRoof;
ParThermalGround_Out(ittm)		=	ParThermalGround;
ParThermalWall_Out(ittm)		=	ParThermalWall;
ParThermalTree_Out(ittm)		=	ParThermalTree;
ParVegRoof_Out(ittm)			=	ParVegRoof;
ParVegGround_Out(ittm)			=	ParVegGround;
ParVegTree_Out(ittm)			=	ParVegTree;

end

Computational_Time =toc;
%profile off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('COMPUTATIONAL TIME [s] ')
disp(Computational_Time)
disp(' COMPUTATIONAL TIME [ms/cycle] ')
disp(1000*Computational_Time/ittn)


% Because we overwrite the initial temperature at first time step we get a
% wrong temperature storage and hence a wrong energy balance in the first
% time step.
for i=1:size(EBNames,1)
	EB.(cell2mat(EBNames(i)))(1,1,:)	=	0;
end

save(['Calculation',NameOutput],'Solver','TempVec','Humidity','SWRabs','SWRin','SWRout','SWREB','LWRabs','LWRin','LWRout',...
	'LWREB','Hflux','LEflux','Gflux','dStorage','RES','Eflux','Runoff','Runon','Leakage',...
	'Int','dInt_dt','Infiltration','Vwater','dVwater_dt','Owater','OSwater','ExWater','SoilPotW',...
	'CiCO2Leaf','WBRoof','WBCanyonIndv','WBCanyonTot','EB','Wind','TempDamp','Qinlat','Results2m',...
	'n','m','Name_Site','MeteoDataRaw','Anthropo',...
	'Gemeotry_m_Out','ParTree_Out','geometry_Out','FractionsRoof_Out','FractionsGround_Out',...
	'WallLayers_Out','ParSoilRoof_Out','ParSoilGround_Out','ParInterceptionTree_Out',...
	'PropOpticalRoof_Out','PropOpticalGround_Out','PropOpticalWall_Out','PropOpticalTree_Out',...
	'ParThermalRoof_Out','ParThermalGround_Out','ParThermalWall_Out','ParThermalTree_Out',...
	'ParVegRoof_Out','ParVegGround_Out','ParVegTree_Out')
