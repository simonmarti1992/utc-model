function[dcan,zomcan,u_Hcan,u_Zp,w_Zp,alpha]=...
	WindProfile_Canyon(Hcan,Htree,R_tree,Wcan,Wroof,Kopt,LAI_t,Zatm,uatm,Zp,trees,Zref_und,zom_und)

% OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dcan		=	Urban displacement height including trees [m]
% zomcan	=	Urban momentum roughness height including trees [m]
% u_Hcan	=	Wind speed at canyon height [m/s]
% u_Zpcan	=	Wind speed within canyon at height Zpcan [m/s]
% w_Zpcan	=	Vertical wind speed within canyon [m/s]

% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hcan		=	canyon height [m]
% Htree		=	Tree height [m]
% R_tree	=	Tree radius [m]
% Wcan		=	Canyon width [m]
% Wroof		=	Roof width [m]
% Kopt		=	Optical transmission factor [-]
% LAI_t		=	Leaf area index of tree [-]	
% Zatm		=	Atmospheric reference height [m]
% uatm		=	Wind speed at atmospheric reference height [m/s]
% Zpcan		=	Height of interest within the canyon [m]
% trees		=	Presence of trees [0: No, 1: Yes]

if trees==0
	Htree	=	0;
	R_tree	=	0;
	LAI_t	=	0;
end

% Displacement height and roughness length according to Kent et al 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters McDonald 1998
a		=	4.43;   % Best fit for staggered arrays , or a = 3.59 for square arrays
k		=	0.4;	% Von Karman costant
b		=	1.0;	% b=1, no incorporation for drag correction factors. Good fit for staggered arrays
CDb		=	1.2;	% nominal drag for cubical obstacles

% Plan area fraction of buildings and vegetation
Ap_build	=	Wroof;
Ap_tree		=	4*R_tree;
Ap_urb		=	Wcan+Wroof;

% Frontal area fraction of vegetation and buildings: assumtpion infinite
% urban canyon perpendicular to the wind direction (Length of building and
% plot equals infinity)
Af_build_s	=	Hcan;
Af_veg_s	=	2*R_tree;

% Tree canopy transmittance (optical = P2D)
P2D			=	exp(-Kopt*LAI_t);
P3D			=	P2D^0.40;			% Guan et al. 2003
Pv			=	(-1.251*P3D^2+ + 0.489*P3D + 0.803)/CDb;	% Guan et al. 2000

% Plan area fraction of buildings and 
% Calculation of structural parameters and wind profile in the city
Lp_tot		=	(Ap_build + (1-P3D)*Ap_tree)/Ap_urb;
H_tot		=	(Hcan*Ap_build + (Htree+R_tree)*(1-P3D)*Ap_tree)/(Ap_build + (1-P3D)*Ap_tree);

% Urban discplacement height and roughness length with incorportation of trees (Kent 2017), (MacDonald
% 1998)
dcan	=	H_tot*(1+a^(-Lp_tot)*(Lp_tot-1));	% displacement height of canyon [m], eq. 23

Af_build=	H_tot/(H_tot-dcan)*Af_build_s;
Af_veg	=	H_tot/(H_tot-dcan)*Af_veg_s;

zomcan	=	H_tot*(1-dcan/H_tot)*exp(-(1/k^2*0.5*b*CDb*(1-dcan/H_tot)*(Af_build+Pv*Af_veg)/Ap_urb)^(-0.5));
zohcan	=	zomcan/10;

% Calculation of wind profile above and in the canyon with a logarithmic
% and exponential wind profile.
us_atm		=	k.*uatm./log((Zatm-dcan)./zomcan);		% Friction Velocity Atmosphere [m/s]
u_Hcan		=	(us_atm./k).*log((Hcan-dcan)./zomcan);	% Wind Speed at canyon top [m/s]
alpha		=	log(uatm./u_Hcan)./(Zatm./Hcan -1);		% Attenuation Coefficient Canyon not corrected for presence of trees.
if Zp >= Hcan
	u_Zp	=	(us_atm./k).*log((Zp-dcan)./zomcan);
	w_Zp	=	0;	
elseif Zp <= Hcan && Zp >= Zref_und
	u_Zp	=	u_Hcan.*exp(-alpha.*(1-Zp./Hcan));
	w_Zp	=	0;	
elseif Zp <= Zref_und && Zp >= zom_und
	uref_und	=	u_Hcan.*exp(-alpha.*(1-Zref_und./Hcan));
	usref_und	=	k*uref_und/log(Zref_und/zom_und);
	u_Zp		=	(usref_und./k).*log(Zp./zom_und);
	w_Zp		=	0;	
else
	u_Zp	=	0;
	w_Zp	=	0;	
	disp('wind speed calculation height higher than reference height or lower than roughness length')
end




