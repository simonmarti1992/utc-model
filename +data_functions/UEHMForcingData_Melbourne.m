function[SunPosition,MeteoData,HumidityAtm,Anthropogenic,location,ParCalculation]...
	=UEHMForcingData(MeteoDataRaw,itt)

% Input variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LWR_in			=	MeteoDataRaw.LWR_in;		%[W/m2]
SAB1_in			=	MeteoDataRaw.SAB1_in;		%[W/m2]
SAB2_in			=	MeteoDataRaw.SAB2_in;		%[W/m2]
SAD1_in			=	MeteoDataRaw.SAD1_in;		%[W/m2]
SAD2_in			=	MeteoDataRaw.SAD2_in;		%[W/m2]
T_atm			=	MeteoDataRaw.T_atm;			%[K]
windspeed_u		=	MeteoDataRaw.windspeed_u;	%[m/s]
pressure_atm	=	MeteoDataRaw.pressure_atm;	%[Pa]
rain			=	MeteoDataRaw.rain;			%[mm/h]
rel_humidity	=	MeteoDataRaw.rel_humidity;	%[-]
date_time		=	MeteoDataRaw.Date;

%% Location properties of the urban area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi				=	dms2degrees([-37,49,0]);	% phi = latitude positive north (degrees)
lambda			=	dms2degrees([144,53,0]);	% lambda = longitude positive west (degrees)
theta_canyon	=	deg2rad(98); % deg2rad(189) % theta_canyon = canyon orientation (rad), measured from google maps: 2 main street directions
DeltaGMT		=	10;							% DeltaGMT = difference with Greenwich Meridian Time [h](local apparent time = UTC+8)(Roth et al. 2016)
		
location		=	struct('phi',phi,'lambda',lambda,'theta_canyon',theta_canyon,...
					'DeltaGMT',DeltaGMT);

%% Calculate zenith angle
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Datam			=	datevec(date_time);
Datam			=	Datam(itt,:);

t_bef			=	0.5;
t_aft			=	0.5;

[h_S,~,zeta_S,~,~,~,~] = data_functions.SetSunVariables(Datam, DeltaGMT, lambda, phi, t_bef,t_aft);

theta_Z			=	pi/2-h_S;	% solar zenith angle

if theta_Z<= -pi/2 || theta_Z>= pi/2
    theta_Z		=	pi/2;
end

theta_n			=	zeta_S-theta_canyon;	% difference between solar azimuth angle and canyon orientation

SunPosition		=	struct('Datam',Datam,'t_bef',t_bef,'t_aft',t_aft,...
					'theta_Z',theta_Z,'theta_n',theta_n);

%% Radiation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SW_dir		=	SAB1_in(itt,1) + SAB2_in(itt,1);	% Direct incoming shortwave radiation W per m^2 of horizontal surface
SW_diff		=	SAD1_in(itt,1) + SAD2_in(itt,1);	% Diffuse incoming shortwave radiation W per m^2 of horizontal surface
LWR			=	LWR_in(itt,1);						% Atmospheric longwave radiation W per m^2 of horizontal surface

if abs(cos(theta_Z))< 0.1	% According to Simones suggestion
    SW_dir		=	0;
end

%% Meteorological data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Zatm			=	40;						% Atmospheric reference height [m]
Tatm			=	T_atm(itt,1);			% Air Temperature at atmospheric reference level [K]
Uatm			=	windspeed_u(itt,1);		% Wind speed at atmospheric reference level [m/s]
Uatm(Uatm==0)	=	0.01;
esat_Tatm		=	611*exp(17.27*(Tatm-273.16)/(237.3+(Tatm-273.16)));	% vapor pressure saturation at Tatm [Pa]
rel_hum			=	rel_humidity(itt,1);	% Relative humidity [-]
ea				=	esat_Tatm*rel_hum;		% vapor pressure [Pa]
Pre				=	pressure_atm(itt,1);	% air pressure [Pa]. Carefull, used to be [hPa - mbar] in T&C
q_atm			=	0.622*ea/(Pre-0.378*ea);% Specifc humidity of air at reference height []
qSat_atm		=	0.622*esat_Tatm/(Pre-0.378*esat_Tatm); % WINDSPEED CANNOT BE 0 OTHERWISE THE LEAF BOUNDARY RESISTANCE FAILS
Catm_CO2		=	400;					% [ppm]-[umolCO2/mol] Atmospheric CO2 concentration 2017
Catm_O2			=	210000;					% [ppm] - [umolO2/mol] Intercellular Partial Pressure Oxygen
Rain			=	rain(itt,1);			% Precipiation [mm]

MeteoData		=	struct('SW_dir',SW_dir,'SW_diff',SW_diff,'LWR',LWR,'Zatm',Zatm,...
					'Tatm',Tatm,'Uatm',Uatm,'esat_Tatm',esat_Tatm,'rel_hum',rel_hum,...
					'ea',ea,'Pre',Pre,'q_atm',q_atm,'Catm_CO2',Catm_CO2,'Catm_O2',Catm_O2,...
					'Rain',Rain);

HumidityAtm		=	struct('AtmRelative',rel_hum,'AtmSpecific',q_atm,'AtmVapourPre',ea,...
					'AtmRelativeSat',1,'AtmSpecificSat',qSat_atm,'AtmVapourPreSat',esat_Tatm);

				

%% ANTHROPOGENIC FACTORS
% Building intertior temperature & anthropogenic heat input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (Tatm-273.15)<18
	Tb = 18+273.15;
elseif (Tatm-273.15)>27
	Tb = 27+273.15;
else	
	Tb	=	Tatm;
end

Qf_canyon	=	0;	% 
Qf_roof		=	0;	% Is so far not considered in the model
% Anthropogenic water
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dth		=	1;	
date				=	datetime(Datam(1),Datam(2),Datam(3),Datam(4),0,0);
Waterf_canyonVeg	=	(date>='15-Nov-2003 00:00:00'& date<'01-Mar-2004 00:00:00')*3/24*dth +...
						(date>='1-Mar-2004 00:00:00'& date<'16-Apr-2004 00:00:00')*2/24*dth;% [mm/h] applied on the vegetated ground surface area
Waterf_canyonBare	=	0;% [mm/h] applied on the vegetated ground surface area
Waterf_roof			=	0;% [mm/h] applied on the vegetated ground surface area

Anthropogenic	=	struct('Tb',Tb,'Qf_canyon',Qf_canyon,'Qf_roof',Qf_roof,...
					'Waterf_canyonVeg',Waterf_canyonVeg,'Waterf_canyonBare',Waterf_canyonBare,'Waterf_roof',Waterf_roof);

%% GENERAL PARAMETERS FOR CALCULATION: DO NOT NEED TO BE CHANGED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dth		=	1;												% time step of calculation [h]
dts		=	3600;											% time step of calculation [s]
row		=	1000;											% density of water [kg/m^3]
cp_atm	=	1005+(((Tatm-273.15)+23.15)^2)/3364;			% specific heat air  [J/kg K]
rho_atm	=	(Pre/(287.04*Tatm))*(1-(ea/Pre)*(1-0.622));		% dry air density at atmosphere [kg/m^3]
					
ParCalculation	=	struct('dts',dts,'dth',dth,'row',row,'cp_atm',cp_atm,'rho_atm',rho_atm);

return


