function[SWRabsRoofImp,SWRabsRoofVeg,SWRabsTotalRoof,...
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
		QRoofImp,QRoofVegDrip,QRoofVegPond,LkRoofImp,LkRoofVeg,LkRoof,QRoofVegSoil,RunoffRoofTot,RunonRoofTot,...
		IntRoofImp,IntRoofVegPlant,IntRoofVegGround,dInt_dtRoofImp,dInt_dtRoofVegPlant,dInt_dtRoofVegGround,...
		IntRooftot,dInt_dtRooftot,dVRoofSoil_dt,...
		fRoofVeg,VRoofSoil,OwRoofSoil,OSwRoofSoil,ExWaterRoof_H,SoilPotWRoof_H,SoilPotWRoof_L,ExWaterRoof_L,...
		CiCO2LeafRoofVegSun,CiCO2LeafRoofVegShd,...
		WBRoofVegInVeg,WBRoofVegInGround,WBRoofVegSoil,...
		EBRoofImp,EBRoofVeg,Yroof,WBRoofImp,WBRoofVeg,WBRoofTot]...
		=EB_WB_roof(TemperatureR,TempVec,MeteoData,itt,Int,ExWater,Vwater,Owater,SoilPotW,CiCO2Leaf,Runon,...
		Gemeotry_m,FractionsRoof,ParSoilRoof,PropOpticalRoof,ParThermalRoof,ParVegRoof,...
		HumidityAtm,Anthropogenic,ParCalculation)

% TemperatureR(1,1) = Troof_imp
% TemperatureR(1,2) = Troof_veg
% TemperatureR(1,3) = Troof_interior_imp
% TemperatureR(1,4) = Troof_interior_veg

%% load parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Gemeotry_m,~,~,FractionsRoof,~,~,ParSoilRoof,~,~,...
% 	PropOpticalRoof,~,~,~,ParThermalRoof,~,~,~,ParVegRoof,~,~]=data_functions.Data_UEHM_site(ittm);
% 
% [~,MeteoData,HumidityAtm,Anthropogenic,~,ParCalculation]...
% 	=data_functions.UEHMForcingData(MeteoDataRaw,itt);

%% Shortwave radiation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SWRabs_dir_veg		=	(1-PropOpticalRoof.aveg)*MeteoData.SW_dir;	% Absorbed direct shortwave radiation by vegetated roof [W/m^2]
SWRabs_diff_veg		=	(1-PropOpticalRoof.aveg)*MeteoData.SW_diff;	% Absorbed diffuse shortwave radiation by vegetated roof [W/m^2]
SWRabsRoofVeg		=	(1-PropOpticalRoof.aveg)*(MeteoData.SW_dir+MeteoData.SW_diff);	% Absorbed total shortwave radiation by vegetated roof [W/m^2]
SWRabsRoofImp		=	(1-PropOpticalRoof.aimp)*(MeteoData.SW_dir+MeteoData.SW_diff);	% Absorbed total shortwave radiation by impervious roof [W/m^2]
SWRabsTotalRoof		=	SWRabsRoofImp*FractionsRoof.fimp + SWRabsRoofVeg*FractionsRoof.fveg;

SWRoutRoofImp		=	PropOpticalRoof.aimp*(MeteoData.SW_dir+MeteoData.SW_diff);
SWRoutRoofVeg		=	PropOpticalRoof.aveg*(MeteoData.SW_dir+MeteoData.SW_diff);
SWRoutTotalRoof		=	SWRoutRoofImp*FractionsRoof.fimp + SWRoutRoofVeg*FractionsRoof.fveg;

SWRinRoofImp		=	(MeteoData.SW_dir+MeteoData.SW_diff);
SWRinRoofVeg		=	(MeteoData.SW_dir+MeteoData.SW_diff);
SWRinTotalRoof		=	(MeteoData.SW_dir+MeteoData.SW_diff);

SWREBRoofImp		=	SWRinRoofImp - SWRoutRoofImp - SWRabsRoofImp;
SWREBRoofVeg		=	SWRinRoofVeg - SWRoutRoofVeg - SWRabsRoofVeg;
SWREBTotalRoof		=	SWRinTotalRoof - SWRoutTotalRoof - SWRabsTotalRoof;

%% Longwave radiation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bolzm				=	5.67*10^(-8);	% Stefan-Boltzmann constant [W*m^-2*K-4], Lecture notes hydrology 2
LWRabsRoofVeg		=	MeteoData.LWR-(PropOpticalRoof.eveg*bolzm*(TemperatureR(1,2))^4+(1-PropOpticalRoof.eveg)*MeteoData.LWR);	% Total absorbed longwave radiation by vegetated roof [W/m^2]
LWRabsRoofImp		=	MeteoData.LWR-(PropOpticalRoof.eimp*bolzm*(TemperatureR(1,1))^4+(1-PropOpticalRoof.eimp)*MeteoData.LWR);	% Total absorbed longwave radiation by impervious roof [W/m^2]
LWRabsTotalRoof		=	LWRabsRoofImp*FractionsRoof.fimp + LWRabsRoofVeg*FractionsRoof.fveg;

LWRoutRoofVeg		=	PropOpticalRoof.eveg*bolzm*(TemperatureR(1,2))^4+(1-PropOpticalRoof.eveg)*MeteoData.LWR;
LWRoutRoofImp		=	PropOpticalRoof.eimp*bolzm*(TemperatureR(1,1))^4+(1-PropOpticalRoof.eimp)*MeteoData.LWR;
LWRoutTotalRoof		=	LWRoutRoofImp*FractionsRoof.fimp + LWRoutRoofVeg*FractionsRoof.fveg;

LWRinRoofImp		=	MeteoData.LWR;
LWRinRoofVeg		=	MeteoData.LWR;
LWRinTotalRoof		=	MeteoData.LWR;

LWREBRoofImp		=	LWRinRoofImp - LWRoutRoofImp - LWRabsRoofImp;
LWREBRoofVeg		=	LWRinRoofVeg - LWRoutRoofVeg - LWRabsRoofVeg;
LWREBTotalRoof		=	LWRinTotalRoof - LWRoutTotalRoof - LWRabsTotalRoof;
	
%% Sensible and latent heat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[HfluxRoofImp,HfluxRoofVeg,EfluxRoofImp,EfluxRoofVegInt,EfluxRoofVegPond,EfluxRoofVegSoil,TEfluxRoofVeg,...
	LEfluxRoofImp,LEfluxRoofVegInt,LEfluxRoofVegPond,LEfluxRoofVegSoil,LTEfluxRoofVeg,...
	CiCO2LeafRoofVegSun,CiCO2LeafRoofVegShd,raRooftoAtm,rb_LRoof,rap_LRoof,r_soilRoof,rs_sunRoof,rs_shdRoof]...
	=turbulent_heat_function.HeatFlux_roof(TemperatureR,MeteoData,HumidityAtm,ParVegRoof,FractionsRoof,Gemeotry_m,...
	ParSoilRoof,ParCalculation,SoilPotW,Owater,Vwater,ExWater,Int,CiCO2Leaf,...
	itt,SWRabs_dir_veg,SWRabs_diff_veg);

HfluxRoof		=	HfluxRoofImp*FractionsRoof.fimp + HfluxRoofVeg*FractionsRoof.fveg;
LEfluxRoofVeg	=	LEfluxRoofVegInt + LEfluxRoofVegPond + LEfluxRoofVegSoil + LTEfluxRoofVeg;
LEfluxRoof		=	LEfluxRoofImp*FractionsRoof.fimp + LEfluxRoofVeg*FractionsRoof.fveg;
EfluxRoofVeg	=	EfluxRoofVegInt + EfluxRoofVegPond + EfluxRoofVegSoil + TEfluxRoofVeg;
EfluxRoof		=	EfluxRoofImp*FractionsRoof.fimp + EfluxRoofVeg*FractionsRoof.fveg;

%% Ground heat flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Impervious conductive heat flux
[G1RoofImp,G2RoofImp,dsRoofImp]=conductive_heat_functions.Impervious_Conductive_HeatRoof...
	(TemperatureR,TempVec,Anthropogenic,ParThermalRoof,ParSoilRoof,ParCalculation,itt);

% Vegetated ground heat flux
[G1RoofVeg,G2RoofVeg,dsRoofVeg]=conductive_heat_functions.Soil_Conductive_Heat...
	(TemperatureR,TempVec,Anthropogenic,Owater,ParVegRoof,ParSoilRoof,ParThermalRoof,ParCalculation,itt);

G1Roof		=	G1RoofImp*FractionsRoof.fimp + G1RoofVeg*FractionsRoof.fveg;
G2Roof		=	G2RoofImp*FractionsRoof.fimp + G2RoofVeg*FractionsRoof.fveg;
dsRoof		=	dsRoofImp*FractionsRoof.fimp + dsRoofVeg*FractionsRoof.fveg;

%% Energy balance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Yroof(1)	=	SWRabsRoofImp+LWRabsRoofImp-HfluxRoofImp-G1RoofImp-LEfluxRoofImp;  % Energy budget impervious roof
Yroof(2)	=	SWRabsRoofVeg+LWRabsRoofVeg-HfluxRoofVeg-G1RoofVeg-LEfluxRoofVegInt-LEfluxRoofVegPond-LEfluxRoofVegSoil-LTEfluxRoofVeg;  % Energy budget vegetated roof
Yroof(3)	=	G1RoofImp-G2RoofImp-dsRoofImp; % Energy budget concrete mass roof
Yroof(4)	=	G1RoofVeg-G2RoofVeg-dsRoofVeg; % Energy budget concrete mass roof

EBRoofImp	=	SWRabsRoofImp+LWRabsRoofImp-HfluxRoofImp-G2RoofImp-dsRoofImp-LEfluxRoofImp;
EBRoofVeg	=	SWRabsRoofVeg+LWRabsRoofVeg-HfluxRoofVeg-G2RoofVeg-dsRoofVeg-LEfluxRoofVegInt-LEfluxRoofVegPond-LEfluxRoofVegSoil-LTEfluxRoofVeg;

%% Water fluxes and water balance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[QRoofImp,IntRoofImp,dInt_dtRoofImp,LkRoofImp,...
	QRoofVegDrip,IntRoofVegPlant,dInt_dtRoofVegPlant,...
	QRoofVegPond,IntRoofVegGround,dInt_dtRoofVegGround,dVRoofSoil_dt,fRoofVeg,...
	VRoofSoil,OwRoofSoil,OSwRoofSoil,LkRoofVeg,SoilPotWRoof_L,ExWaterRoof_L,QRoofVegSoil,TEfluxRoofVeg,EfluxRoofVegSoil,RunoffRoofTot,RunonRoofTot,...
	~,WBRoofVegInVeg,WBRoofVegInGround,WBRoofVegSoil,...
	WBRoofImp,WBRoofVeg,WBRoofTot]=...
	water_functions.WaterRoof(EfluxRoofImp,EfluxRoofVegInt,EfluxRoofVegPond,EfluxRoofVegSoil,TEfluxRoofVeg,...
	MeteoData,Int,Owater,Runon,FractionsRoof,ParSoilRoof,ParCalculation,ParVegRoof,Anthropogenic,itt);

LkRoof			=	LkRoofImp*FractionsRoof.fimp + LkRoofVeg*FractionsRoof.fveg;
IntRooftot		=	IntRoofImp*FractionsRoof.fimp + (IntRoofVegPlant+IntRoofVegGround)*FractionsRoof.fveg;
dInt_dtRooftot	=	dInt_dtRoofImp*FractionsRoof.fimp + (dInt_dtRoofVegPlant+dInt_dtRoofVegGround)*FractionsRoof.fveg;
ExWaterRoof_H	=	NaN(1,ParSoilRoof.ms);
SoilPotWRoof_H	=	NaN(1,1);

end


