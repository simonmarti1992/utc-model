function[Gemeotry_m,ParTree,geometry,FractionsRoof,FractionsGround,...
	WallLayers,ParSoilRoof,ParSoilGround,ParInterceptionTree,...
	PropOpticalRoof,PropOpticalGround,PropOpticalWall,PropOpticalTree,...
	ParThermalRoof,ParThermalGround,ParThermalWall,ParThermalTree,...
	ParVegRoof,ParVegGround,ParVegTree]=Data_UEHM_site(varargin)

%% GEOMETRY OF URBAN AREA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Height_canyon	=	4.5;	% Canyon (building) height (m)
Width_canyon	=	11.3;	% Canyon (road) width (m)
Width_roof		=	4;		% Roof width (m)
Height_tree		=	4;		% Tree height (m)
Radius_tree		=	0.19;	% Tree-crown radius (m)
Distance_tree	=	2;		% Tree-to-wall distance (m)

trees	=	1;		% If trees are present in the urban canyon: 1, if no trees are present: 0
ftree	=	1;		% Tree fraction along canyon axis (Always keep as 1)

hcanyon			=	Height_canyon/Width_canyon;		% normalized canyon height(-)
wcanyon			=	Width_canyon/Width_canyon;		% normalized canyon width (-)
wroof			=	Width_roof/Width_canyon;		% normalized roof width (-)
htree			=	Height_tree/Width_canyon;		% normalized tree height (-)
radius_tree		=	Radius_tree/Width_canyon;		% normalized tree radius (-)
distance_tree	=	Distance_tree/Width_canyon;		% normalized tree-to-wall distance (-)
ratio			=	hcanyon/wcanyon;				% height to width ratio (-)

wcanyon_norm	=	wcanyon/(wcanyon+wroof);		% normalized canyon width overall (-)
wroof_norm		=	wroof/(wcanyon+wroof);			% normalized roof width overall (-)

Gemeotry_m	=	struct('Height_canyon',Height_canyon,'Width_canyon',Width_canyon,...
				'Width_roof',Width_roof,'Height_tree',Height_tree,...
				'Radius_tree',Radius_tree,'Distance_tree',Distance_tree);

ParTree		=	struct('trees',trees,'ftree',ftree);

geometry	=	struct('hcanyon',hcanyon,'wcanyon',wcanyon,'wroof',wroof,...
				'htree',htree,'radius_tree',radius_tree,'distance_tree',distance_tree,...
				'ratio',ratio,'wcanyon_norm',wcanyon_norm,'wroof_norm',wroof_norm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SURFACE FRACTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fveg_R	=	0;	% Partitioning roof vegetation surface(-)
fimp_R	=	1;	% Partitioning roof impervious(-)

if fveg_R+fimp_R~=1
	disp('roof fractions do not add up to 1, check Data_UEHM_roof.m')
end

Per_runoff_R	=	1;		% Percentage of excess water that leaves the roof system as runoff [-]

FractionsRoof	=	struct('fveg',fveg_R,'fimp',fimp_R,'Per_runoff',Per_runoff_R);

% GROUND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fveg_G	=	0.15;	% Partitioning ground vegetation surface(-)
fbare_G	=	0.53;	% Partitioning ground bare surface(-)
fimp_G	=	0.32;	% Partitioning ground impervious(-)

if fveg_G+fbare_G+fimp_G~=1
	disp('ground fractions do not add up to 1, check Data_UEHM_ground.m')
end

Per_runoff_G	=	0.5;	% Percentage of excess water that leaves the canyon system as runoff [-]

FractionsGround	=	struct('fveg',fveg_G,'fbare',fbare_G,'fimp',fimp_G,'Per_runoff',Per_runoff_G);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTICAL PROPERTIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROOF OPTICAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aveg_R		=	0.27;	% Roof vegetation surface albedo (-)
aimp_R		=	0.16;	% Roof impervious albedo (-)
albedo_R	=	fveg_R*aveg_R+fimp_R*aimp_R;	% Equivalent roof surface albedo (-)

eveg_R		=	0.97;	% Roof vegetation surface emissivity (-)
eimp_R		=	0.95;	% Roof impervious emissivity (-)
emissivity_R=	fveg_R*eveg_R+fimp_R*eimp_R;	% Equivalent roof surface emissivity (-)

PropOpticalRoof	=	struct('aveg',aveg_R,'aimp',aimp_R,'albedo',albedo_R,...
					'eveg',eveg_R,'eimp',eimp_R,'emissivity',emissivity_R);

% GROUND OPTICAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aveg_G		=	0.27;	% Ground vegetation surface albedo (-)
abare_G		=	0.2;	% Ground bare surface albedo (-)
aimp_G		=	0.15;	% Ground impervious albedo (-)
albedo_G	=	fveg_G*aveg_G + fbare_G*abare_G + fimp_G*aimp_G;	% equivalent ground surface albedo (-)

eveg_G		=	0.97;	% Ground vegetation surface emissivity (-)
ebare_G		=	0.98;	% Ground bare surface emissivity (-)
eimp_G		=	0.95;	% Ground impervious emissivity (-)
emissivity_G=	fveg_G*eveg_G + fbare_G*ebare_G + fimp_G*eimp_G;	% equivalent ground surface emissivity (-)

PropOpticalGround	=	struct('aveg',aveg_G,'abare',abare_G,'aimp',aimp_G,'albedo',albedo_G,...
						'eveg',eveg_G,'ebare',ebare_G,'eimp',eimp_G,'emissivity',emissivity_G);
					
% WALL OPTICAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
albedo_W		=	0.5;	% Wall surface albedo (-)
emissivity_W	=	0.95;	% Wall emissivity (-)

PropOpticalWall	=	struct('albedo',albedo_W,'emissivity',emissivity_W);

% TREE OPTICAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
albedo_T		=	0.27;	% Tree albedo (-)
emissivity_T	=	0.97;	% Tree emissivity (-)

PropOpticalTree	=	struct('albedo',albedo_T,'emissivity',emissivity_T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% THERMAL PROPERTIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROOF THERMAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lan_dry_imp_R	=	0.6;		% Thermal conductivity dry solid, impervious roof [W/m K]
cv_s_imp_R		=	1.9*10^6;	% Volumetric heat capacity solid, impervious roof [J/m^3 K]

ParThermalRoof	=	struct('lan_dry_imp',lan_dry_imp_R,'cv_s_imp',cv_s_imp_R);

% GROUND THERMAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lan_dry_imp_G	=	1.2;		% Thermal conductivity dry solid, impervious ground [W/m K]
cv_s_imp_G		=	1.1*10^6;	% Volumetric heat capacity solid, impervious ground [J/m^3 K]

ParThermalGround	=	struct('lan_dry_imp',lan_dry_imp_G,'cv_s_imp',cv_s_imp_G);

% WALL THERMAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lan_dry_W	=	1.3;			% Thermal conductivity dry solid, wall [W/m K]
cv_s_W		=	1.5*10^6;		% Volumetric heat capacity solid, wall [J/m^3 K]

ParThermalWall	=	struct('lan_dry',lan_dry_W,'cv_s',cv_s_W);

% VEGETATION THERMAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cthermal_leaf	=	640;		% [J m-2 K-1] Heat capacity per single leaf area based on Kitaya et al. 2003, Ryu et al. 2016
								% Not used in this version

ParThermalTree	=	struct('Cthermal_leaf',Cthermal_leaf);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROOF VEGETATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General
LAI_R		=	4;			% Leaf area index, roof vegetation (-)
SAI_R		=	0;			% Stem area index, roof vegetation (-)
hc_R		=	0.5;		% canopy height, roof vegetation
h_disp_R	=	2/3*hc_R;	% Zero plane displacement height, roof vegetation [m]
d_leaf_R	=	5;			% Leaf dimension, roof vegetation [cm]

% Roof water uptake
CASE_ROOT_R	=	1;		% Type of Root Profile, roof vegetation [-]
ZR95_R		=	1;		% Root depth 95 percentile, roof vegetation [mm]
ZR50_R		=	NaN;	% Root depth 50 percentile, roof vegetation [mm]
ZRmax_R		=	NaN;	% Maximum Root depth, roof vegetation [mm]
Rrootl_R	=	3000;	% Root length index, roof vegetation [m root / m^2 PFT]
PsiL50_R	=	-4.0;	% Water Potential at 50% loss conductivity, roof vegetation [MPa]  
PsiX50_R	=	-4.5;	% Water potential at 50 of xylem hydraulic conductivity and limit for water extraction from soil, roof vegetation [MPa]

% Photosynthesis and Transpiration
FI_R		=	0.081;	% Intrinsec quantum Efficiency, roof vegetation [umolCO2/umolPhotons]
Do_R		=	1000;	% Empirical coefficient for the role of vapor pressure in the biochemical model of photosynthesis, roof vegetation [Pa] 
a1_R		=	6;		% Empirical parameter connecting stomatal aperture and net assimilaton, roof vegetation [-]
go_R		=	0.01;	% minimal Stomatal Conductance,roof vegetation [mol / s m^2] 
CT_R		=	3;		% Photosyntesis Typology for Plants, Photosynthetic pathway C3 or C4 --> 'CT' == 3  'CT' ==  4, roof vegetation
DSE_R		=	0.656;	% Activation Energy in Photosynthesis for Rubisco Capacity, roof vegetation [kJ/mol] 
Ha_R		=	55;		% Activation energy, roof vegetation [kJ / mol K]  
gmes_R		=	Inf;	% Mesophyll conductance, not used [mol CO2/s m2] 
rjv_R		=	2.4;	% Scaling factor between Jmax and Vmax, roof vegetation [mol Eq/ mol CO2]
Kopt_R		=	0.5;	% optical depth of direct beam perunit plant area, roof vegetation [-]
Knit_R		=	0.15;	% Canopy nitrogen decay coefficient, roof vegetation [-]
Vmax_R		=	96;		% Maximum Rubisco capacity at 25°C leaf level, roof vegetation [mol CO2/ m2 s]
mSl_R		=	0.0;	% Not used
e_rel_R		=	1;		% Relative Efficiency of the photosynthesis apparatus due to Age/Day-length, roof vegetation [-]
e_relN_R	=	1;		% Relative efficiency of the photosynthesis apparatus due to N limitations, roof vegetation [-]
Psi_sto_00_R=	-0.5;	% Water Potential at PLCs loss conductivity, roof vegetation [MPa]
Psi_sto_50_R=	-3.0;	% Water Potential at 50% loss conductivity, roof vegetation [MPa]
Sl_R		=	0.035;	% Specific leaf area of  biomass, roof vegetation [m^2 /gC]

ParVegRoof	=	struct('LAI',LAI_R,'SAI',SAI_R,'hc',hc_R,'h_disp',h_disp_R,...
				'd_leaf',d_leaf_R,'CASE_ROOT',CASE_ROOT_R,'ZR95',ZR95_R,'ZR50',ZR50_R,...
				'ZRmax',ZRmax_R,'Rrootl',Rrootl_R,'PsiL50',PsiL50_R,'PsiX50',PsiX50_R,...
				'FI',FI_R,'Do',Do_R,'a1',a1_R,'go',go_R,'CT',CT_R,'DSE',DSE_R,'Ha',Ha_R,...
				'gmes',gmes_R,'rjv',rjv_R,'Kopt',Kopt_R,'Knit',Knit_R,'Vmax',Vmax_R,...
				'mSl',mSl_R,'e_rel',e_rel_R,'e_relN',e_relN_R,'Psi_sto_00',Psi_sto_00_R,...
				'Psi_sto_50',Psi_sto_50_R,'Sl',Sl_R);
				
% GROUND VEGETATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRASS
% General
LAI_G		=	1.5;		% Leaf area index of ground vegetation (-)
SAI_G		=	0.001;		% Stem area index of ground vegetation (-)
hc_G		=	0.1;		% canopy height of ground vegetation [m]
h_disp_G	=	2/3*hc_G;	% Zero plane displacement height of ground vegetation [m]
d_leaf_G	=	0.8;		% Leaf dimension of ground vegetation [cm]

% Ground water uptake
CASE_ROOT_G	=	1;		% Type of Root Profile, ground vegetation [-]
ZR95_G		=	250;	% Root depth 95 percentile, ground vegetation [mm]
ZR50_G		=	NaN;	% Root depth 50 percentile, ground vegetation [mm]
ZRmax_G		=	NaN;	% Maximum Root depth, ground vegetation [mm]
Rrootl_G	=	2000;	% Root length index, ground vegetation [m root / m^2 PFT]
PsiL50_G	=	-2.5;	% Water Potential at 50% loss conductivity, ground vegetation [MPa]
PsiX50_G	=	-3.5;	% Water potential at 50 of xylem hydraulic conductivity and limit for water extraction from soil, ground vegetation [MPa]

% Photosynthesis and Transpiration
FI_G		=	0.081;	% Intrinsec quantum Efficiency, ground vegetation [umolCO2/umolPhotons]
Do_G		=	2000;	% Empirical coefficient for the role of vapor pressure in the biochemical model of photosynthesis, ground vegetation [Pa] 
a1_G		=	6;		% Empirical parameter connecting stomatal aperture and net assimilaton, ground vegetation [-]
go_G		=	0.01;	% minimal Stomatal Conductance,ground vegetation [mol / s m^2]
CT_G		=	3;		% Photosyntesis Typology for Plants, Photosynthetic pathway C3 or C4 --> 'CT' == 3  'CT' ==  4, ground vegetation
DSE_G		=	0.649;	% Activation Energy in Photosynthesis for Rubisco Capacity, ground vegetation [kJ/mol]
Ha_G		=	72;		% Activation energy, ground vegetation [kJ / mol K]
gmes_G		=	Inf;	% Mesophyll conductance, not used [mol CO2/s m2] 
rjv_G		=	2.2;	% Scaling factor between Jmax and Vmax, ground vegetation [mol Eq/ mol CO2]
Kopt_G		=	0.5;	% optical depth of direct beam perunit plant area, ground vegetation [-]
Knit_G		=	0.2;	% Canopy nitrogen decay coefficient, ground vegetation [-]
Vmax_G		=	58;		% Maximum Rubisco capacity at 25°C leaf level, ground vegetation [mol CO2/ m2 s]
mSl_G		=	0.0;	% Not used
e_rel_G		=	1;		% Relative Efficiency of the photosynthesis apparatus due to Age/Day-length, ground vegetation [-]
e_relN_G	=	1;		% Relative efficiency of the photosynthesis apparatus due to N limitations, ground vegetation [-]
Psi_sto_00_G=	-0.5;	% Water Potential at PLCs loss conductivity, ground vegetation [MPa]
Psi_sto_50_G=	-3;		% Water Potential at 50% loss conductivity, ground vegetation [MPa]
Sl_G		=	0.022;	% Specific leaf area of  biomass, ground vegetation [m^2 /gC]

ParVegGround	=	struct('LAI',LAI_G,'SAI',SAI_G,'hc',hc_G,'h_disp',h_disp_G,...
					'd_leaf',d_leaf_G,'CASE_ROOT',CASE_ROOT_G,'ZR95',ZR95_G,'ZR50',ZR50_G,...
					'ZRmax',ZRmax_G,'Rrootl',Rrootl_G,'PsiL50',PsiL50_G,'PsiX50',PsiX50_G,...
					'FI',FI_G,'Do',Do_G,'a1',a1_G,'go',go_G,'CT',CT_G,'DSE',DSE_G,'Ha',Ha_G,...
					'gmes',gmes_G,'rjv',rjv_G,'Kopt',Kopt_G,'Knit',Knit_G,'Vmax',Vmax_G,...
					'mSl',mSl_G,'e_rel',e_rel_G,'e_relN',e_relN_G,'Psi_sto_00',Psi_sto_00_G,...
					'Psi_sto_50',Psi_sto_50_G,'Sl',Sl_G);
					
% TREE VEGETATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General
LAI_T		=	1.8;		% Leaf area index for tree (-)
SAI_T		=	0.1;		% Stem area index for tree(-)
d_leaf_T	=	1.5;		% Leaf dimension of tree [cm]

% Tree water uptake
CASE_ROOT_T	=	1;		% Type of Root Profile, tree [-]
ZR95_T		=	1000;	% Root depth 95 percentile, tree [mm]
ZR50_T		=	NaN;	% Root depth 50 percentile, tree [mm]
ZRmax_T		=	NaN;	% Maximum Root depth, tree [mm]
Rrootl_T	=	1200;	% Root length index, tree [m root / m^2 PFT]
PsiL50_T	=	-1.2;	% Water Potential at 50% loss conductivity, tree [MPa]
PsiX50_T	=	-4;		% Water potential at 50 of xylem hydraulic conductivity and limit for water extraction from soil, tree [MPa]

% Photosynthesis and Transpiration
FI_T		=	0.081;	% Intrinsec quantum Efficiency, tree [umolCO2/umolPhotons]
Do_T		=	2000;	% Empirical coefficient for the role of vapor pressure in the biochemical model of photosynthesis, tree [Pa] 
a1_T		=	9;		% Empirical parameter connecting stomatal aperture and net assimilaton, tree [-]
go_T		=	0.01;	% minimal Stomatal Conductance, tree [mol / s m^2]
CT_T		=	3;		% Photosyntesis Typology for Plants, Photosynthetic pathway C3 or C4 --> 'CT' == 3  'CT' ==  4, tree
DSE_T		=	0.649;	% Activation Energy in Photosynthesis for Rubisco Capacity, tree [kJ/mol]
Ha_T		=	72;		% Activation energy, tree [kJ / mol K]
gmes_T		=	Inf;	% Mesophyll conductance, not used [mol CO2/s m2] 
rjv_T		=	2;		% Scaling factor between Jmax and Vmax, tree [mol Eq/ mol CO2]
Kopt_T		=	0.5;	% optical depth of direct beam perunit plant area, tree [-]
Knit_T		=	0.25;	% Canopy nitrogen decay coefficient, tree [-]
Vmax_T		=	45;		% Maximum Rubisco capacity at 25°C leaf level, tree [mol CO2/ m2 s]
mSl_T		=	0.0;	% Not used
e_rel_T		=	1;		% Relative Efficiency of the photosynthesis apparatus due to Age/Day-length, tree [-]
e_relN_T	=	1;		% Relative efficiency of the photosynthesis apparatus due to N limitations, tree [-]
Psi_sto_00_T=	-0.9;	% Water Potential at PLCs loss conductivity, tree [MPa]
Psi_sto_50_T=	-2;		% Water Potential at 50% loss conductivity, tree [MPa]
Sl_T		=	0.015;	% Specific leaf area of  biomass, tree [m^2 /gC]
SPARTREE	=	2;		% Tree root distribution: 1 = Tree roots can access all water in the soil (imp, bare, veg) equally
						% 2 =  If the tree crown is smaller than the combined vegetated and bare fraction, 
						% then the trees only transpire from these fractions. Otherwise, they
						% also transpire from the impervious ground fraction.

ParVegTree	=	struct('LAI',LAI_T,'SAI',SAI_T,'d_leaf',d_leaf_T,'CASE_ROOT',CASE_ROOT_T,...
				'ZR95',ZR95_T,'ZR50',ZR50_T,'ZRmax',ZRmax_T,'Rrootl',Rrootl_T,...
				'PsiL50',PsiL50_T,'PsiX50',PsiX50_T,'FI',FI_T,'Do',Do_T,'a1',a1_T,...
				'go',go_T,'CT',CT_T,'DSE',DSE_T,'Ha',Ha_T,'gmes',gmes_T,'rjv',rjv_T,...
				'Kopt',Kopt_T,'Knit',Knit_T,'Vmax',Vmax_T,'mSl',mSl_T,'e_rel',e_rel_T,...
				'e_relN',e_relN_T,'Psi_sto_00',Psi_sto_00_T,'Psi_sto_50',Psi_sto_50_T,...
				'Sl',Sl_T,'SPARTREE',SPARTREE);
					
					
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOLID LAYER DISCRETIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Layer thickness
dz1_R	=	0.075;	% Thickness of first roof layer [m]
dz2_R	=	0.075;	% Thickness of second roof layer [m]

% Soil layer discretization 
Zs_R	=	[ 0 10 20 30 75];	% Roof soil layer discretization [mm]
ms_R	=	length(Zs_R)-1;		% Number of soil layers, roof [-]

% GROUND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Soil layer discretization 
Zs_G	=	[0 10 20 50 100 150 200 300 400 600 800 1000 1500 2000];	% Ground soil layer discretization [mm]
ms_G	=	length(Zs_G)-1;												% number of soil layers, groun [-]

% WALL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Layer thickness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dz1_W	=	0.075;	% Thickness of first wall layer [m]
dz2_W	=	0.075;	% Thickness of second wall layer [m]

WallLayers	=	struct('dz1_wall',dz1_W,'dz2_wall',dz2_W);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INTERCEPTION AND SOIL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Soil classification: Sandy Loam
Pcla_R	=	0.20;    % Fraction of clay in the soil, vegetated roof [-]
Psan_R	=	0.40;    % Fraction of sand in the soil, vegetated roof [-]
Porg_R	=	0.025;   % Fraction of organic material in the soil, vegetated roof [-]

% Interception and soil parameters
In_max_imp_R	=	0.25;	% Maxiumum interception capacity of roof impervious area [mm]
In_max_ground_R	=	10;		% Maxiumum interception capacity of ground under roof vegetation [mm]
Sp_In_R			=	0.2;	% specific water retained by a vegetated roof surface [mm m^2 VEG area m^-2 plant area]

Kimp_R	=	0;		% Hydraulic conductivity of impervious roof area [mm/h]
Kfc_R	=	0.2;	% Conductivity at field capacity, vegetated roof [mm/h]
Phy_R	=	10000;	% Suction at the residual/hygroscopic water content, vegetated roof [kPa]
SPAR_R	=	2;		% SOIL PARAMETER TYPE 1-VanGenuchten 2-Saxton-Rawls, vegetated roof
Kbot_R	=	NaN;	% Conductivity at the bedrock layer, vegetated roof [mm/h] 

ParSoilRoof	=	struct('Zs',Zs_R,'ms',ms_R,'dz1',dz1_R,'dz2',dz2_R,...
				'In_max_imp',In_max_imp_R,'In_max_ground',In_max_ground_R,...
				'Sp_In',Sp_In_R,'Kimp',Kimp_R,'Kfc',Kfc_R,'Phy',Phy_R,'SPAR',SPAR_R,...
				'Kbot',Kbot_R,'Pcla',Pcla_R,'Psan',Psan_R,'Porg',Porg_R);

% GROUND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Soil classification: Sandy Loam
Pcla_G	=	0.20;	% Fraction of clay in the soil [-]
Psan_G	=	0.40;	% Fraction of sand in the soil [-]
Porg_G	=	0.025;	% Fraction of organic material in the soil [-]

% Interception and soil parameters
In_max_imp_G		=	0.5;	% Maxiumum interception capacity of impervious ground area [mm]
In_max_underveg_G	=	10;		% Maxiumum interception capacity of vegetated ground area [mm]
In_max_bare_G		=	10;		% Maxiumum interception capacity of bare ground area [mm]
Sp_In_G				=	0.2;	% specific water retained by a vegetated surface on the ground [mm m^2 VEG area m^-2 plant area]

Kimp_G	=	0.001;		% Hydraulic conductivity of impervious ground area [mm/h]

Kfc_G	=	0.2;		% Conductivity at field capacity of ground [mm/h]
Phy_G	=	10000;		% Suction at the residual/hygroscopic water content of ground [kPa]
SPAR_G	=	2;			% SOIL PARAMETER TYPE 1-VanGenuchten 2-Saxton-Rawls, ground
Kbot_G	=	NaN;		% Conductivity at the bedrock layer, ground [mm/h] 

ParSoilGround	=	struct('Zs',Zs_G,'ms',ms_G,'In_max_imp',In_max_imp_G,...
					'In_max_underveg',In_max_underveg_G,'In_max_bare',In_max_bare_G,...
					'Sp_In',Sp_In_G,'Kimp',Kimp_G,'Kfc',Kfc_G,'Phy',Phy_G,'SPAR',SPAR_G,...
					'Kbot',Kbot_G,'Pcla',Pcla_G,'Psan',Psan_G,'Porg',Porg_G);

% TREES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interception tree
Sp_In_T		=	0.1;	% specific water retained by the tree [mm m^2 VEG area m^-2 plant area]

ParInterceptionTree	=	struct('Sp_In',Sp_In_T);


