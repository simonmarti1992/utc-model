
Hcan		=	30;
Htree		=	[5, 10, 15, 20];
R_tree		=	1.5;
Wcan		=	20;
Wroof		=	10;
Kopt		=	0.61;
LAI_t		=	4;
Zatm		=	40;
uatm		=	10;
Zp			=	[0:0.1:Zatm];
trees		=	1;
Zref_und	=	1.5;
zom_und		=	0.123*0.2;

u_Zp		=	NaN(length(Zp),4);
w_Zp		=	NaN(length(Zp),4);

for j=1:4
	for i=1:length(Zp)
	[dcan,zomcan,u_Hcan,u_Zp(i,j),w_Zp(i,j)]=resistance_functions.WindProfile_Canyon...
		(Hcan,Htree(j),R_tree,Wcan,Wroof,Kopt,LAI_t,Zatm,uatm,Zp(i),trees,Zref_und,zom_und);
	end
end

figure
plot(u_Zp(:,1)./uatm,Zp,'DisplayName','Wind speed above canyon[m/s]')
hold on
plot(u_Zp(:,2)./uatm,Zp,'DisplayName','Wind speed above canyon[m/s]')
plot(u_Zp(:,3)./uatm,Zp,'DisplayName','Wind speed above canyon[m/s]')
plot(u_Zp(:,4)./uatm,Zp,'DisplayName','Wind speed above canyon[m/s]')
ylabel('Canyon height [m]')
xlabel('mean wind speed u [m/s]')
legend('show')
