function[Hwsun,Hwshade,Ewsun,Ewshade,LEwsun,LEwshade,RES_w1,RES_w2,rap_Zp1_In,rap_Zp2_In,...
	Hwsun1,Hwshade1,Hwsun2,Hwshade2,cp_atm,rho_atm,L_heat,Zp1,Zp2,rap_Zp1,rap_Zp2]...
	=HeatFlux_wall(TemperatureC,Gemeotry_m,MeteoData,ParVegTree,ParTree,ParVegGround,FractionsGround)


%% Parameter definitions
Twsun			=	TemperatureC(1,4);
Twshade			=	TemperatureC(1,5);
T_canyon		=	TemperatureC(1,9);
H				=	Gemeotry_m.Height_canyon;
W				=	Gemeotry_m.Width_canyon;
Wroof			=	Gemeotry_m.Width_roof;
Htree			=	Gemeotry_m.Height_tree;
R_tree			=	Gemeotry_m.Radius_tree;
Kopt			=	ParVegTree.Kopt;
LAI_t			=	ParVegTree.LAI;
trees			=	ParTree.trees;
Zatm			=	MeteoData.Zatm;
Tatm			=	MeteoData.Tatm;
Uatm			=	MeteoData.Uatm;
Pre				=	MeteoData.Pre;
ea				=	MeteoData.ea;
hc_L			=	ParVegGround.hc;
fgveg			=	FractionsGround.fveg;
fgbare			=	FractionsGround.fbare;
fgimp			=	FractionsGround.fimp;

%% OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hwsun = sensible heat of sunlit wall (W/m^2)
% Hwshade = sensible heat of shaded wall (W/m^2)
% Ewsun = latent heat of sunlit wall (x)
% Ewshade = latent heat of shaded wall (x)

% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Twsun = Temperature of the wall sun [K]
% Twshade = Temperature of the wall shade [K]
% T_canyon = Temperature of canyon air [K]
% cp_atm = specific heat air  [J/kg K]
% rho_atm = dry air density at atmosphere [kg/m^3]
% U_canyon = horizontal wind speed in canyon at H/2 [m/s]
% W_canyon = Vertical wind speed should be 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensible heat flux sunlit and shaded wall
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cp_atm		=	1005+(((Tatm-273.15)+23.15)^2)/3364;		% specific heat air  [J/kg K]
rho_atm		=	(Pre/(287.04*Tatm))*(1-(ea/Pre)*(1-0.622));	% dry air density at atmosphere [kg/m^3]
% q_atm		=	0.622*ea/(Pre-0.378*ea);					% Specifc humidity of air at reference height []
L_heat		=	1000*(2501.3 - 2.361*(Tatm-273.15));				% Latent heat vaporization/condensaition [J/kg]

% Calcualtion of structural parameters and wind profile in the city
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,~,zom_ground,~,~,~,~,~,~,~,~,~]...
	=resistance_functions.Urban_roughness(0,(fgveg>0)*hc_L,(fgbare>0),(fgimp>0),0);

[dcan,zomcan,~,~,~,~]=resistance_functions.WindProfile_Canyon...
	(H,Htree,R_tree,W,Wroof,Kopt,LAI_t,Zatm,Uatm,2,trees,1.5,zom_ground);

Zp1	=	2;
Zp2	=	2*Zp1+(H-2*Zp1)/2;
wcan=	0;

[~,rap_Zp1,rap_Zp1_In,rap_Zp2,rap_Zp2_In,~,~,~,u_Zp1,u_Zp2,~,~,~]...
	=resistance_functions.InCanyonAerodynamicResistance(Uatm,Zatm,T_canyon-273.15,(Twsun+Twshade)/2-273.15,...
	H,dcan,zomcan,1.5,zom_ground,Zp1,Zp2,2);

RES_w1	=	cp_atm*rho_atm*(11.8+4.2*sqrt((u_Zp1)^2+(wcan)^2))^(-1);
RES_w2	=	cp_atm*rho_atm*(11.8+4.2*sqrt((u_Zp2)^2+(wcan)^2))^(-1);

Hwsun1		=	cp_atm*rho_atm*(Twsun-T_canyon)/(RES_w1+rap_Zp1_In);
Hwshade1	=	cp_atm*rho_atm*(Twshade-T_canyon)/(RES_w1+rap_Zp1_In); 

Hwsun2		=	cp_atm*rho_atm*(Twsun-T_canyon)/(RES_w2+rap_Zp2_In);
Hwshade2	=	cp_atm*rho_atm*(Twshade-T_canyon)/(RES_w2+rap_Zp2_In);

Hwsun		=	max(2*Zp1/H,1)*Hwsun1 + min((H-2*Zp1)/H,0)*Hwsun2;
Hwshade		=	max(2*Zp1/H,1)*Hwshade1 + min((H-2*Zp1)/H,0)*Hwshade2;

Ewsun		=	0;
Ewshade		=	0;
LEwsun		=	L_heat*Ewsun;
LEwshade	=	L_heat*Ewshade; 
