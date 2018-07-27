function plotLH2Data(data)
% plotLH2Data(data)
%	Plots results from data.
%	'data' is a data structure returned from LH2Simulate.


% here : (ST) or 1 is 17,000 gallon horizontal trailer - feeding vessel 
%        (ET) or 2 is 3,300 gallon vertical storage - receiving vessel
%   

try
	P = evalin('base','LH2Model');
catch ME
	if strcmp(ME.identifier,'MATLAB:UndefinedFunction')
		evalin('base','LH2ModelParams');
		P = evalin('base','LH2Model');
	else
		error(ME.message);
	end
end



% liquid levels in (ST) and (ET) and transfer flow, figure 1
figure;
subplot(2,2,1);
plot(data.t/60,data.hL1)
ylabel('(ST) Height (m)');
xlim([0 data.t(end)/60]);
xlabel('Time (min)');
grid on;

subplot(2,2,2);
plot(data.t/60,data.hL2,data.t/60,0.70*P.H*ones(size(data.t)),data.t/60,0.80*P.H*ones(size(data.t)),data.t/60,0.90*P.H*ones(size(data.t)),data.t/60,P.H*ones(size(data.t)))
ylabel('(ET) Height (m)');
xlim([0 data.t(end)/60]);
ylim([0 4.1]);
xlabel('Time (min)');
legend('level','70pct','80pct','90pct','100pct');
grid on;

subplot(2,2,3:4);
plot(data.t/60,data.Jtr);
ylabel('Transfer Line Mass Flow (kg/s)');
xlim([0 data.t(end)/60]);
legend('Jtr');
xlabel('Time (min)');
grid on;

% temperatures,  figure 2
figure;
subplot(2,2,1:2);
plot(data.t/60,data.TL1(:,end),data.t/60,data.Ts1,data.t/60,data.Tv1(:,end))
ylabel('Temperatures in (ST) (K)');
xlim([0 data.t(end)/60]);
xlabel('Time (min)');
legend('Liquid','Surface','Vapor');
grid on;

subplot(2,2,3:4);
plot(data.t/60,data.TL2(:,end),data.t/60,data.Ts2,data.t/60,data.Tv2(:,end),data.t/60,data.Tw2)
ylabel('Temperatures in (ET) (K)');
xlim([0 data.t(end)/60]);
xlabel('Time (min)');
legend('Liquid','Surface','Vapor','Twall');
grid on;

% pressures and internal energies, figure 3
figure;
subplot(2,2,1);
plot(data.t/60,data.uv1(:,end),data.t/60,data.uL1(:,end),data.t/60,data.uv2(:,end),data.t/60,data.uL2(:,end))
ylabel('Specific internal energies (kJ/kg)');
xlim([0 data.t(end)/60]);
xlabel('Time (min)');
legend('uv1','uL1','uv2','uL2');
grid on;

subplot(2,2,2);
plot(data.t/60,data.pv1/6894.75729,data.t/60,data.pv2/6894.75729,data.t/60,(P.p_ET_high/6894.75729)*ones(size(data.t)),data.t/60,((P.p_ET_low/6894.75729))*ones(size(data.t)))
ylabel('Pressure (psia)');
xlabel('Time (min)');
legend('Pv1','Pv2','Upper Threshold (ET)','Lower Threshold (ET)');
grid on;
xlim([0 data.t(end)/60]);

subplot(2,2,3);
plot(data.t/60,data.mL1,data.t/60,data.mL2);
ylabel('mass of liquid (kg)');
xlim([0 data.t(end)/60]);
xlabel('Time (min)');
legend('(ST)','(ET)');
grid on;
xlim([0 data.t(end)/60]);

subplot(2,2,4);
plot(data.t/60,data.mv1,data.t/60,data.mVap,data.t/60,data.mv2);
ylabel('mass of vapor (kg)');
xlim([0 data.t(end)/60]);
xlabel('Time (min)');
legend('(ST)','vaporizer','(ET)');
grid on;
xlim([0 data.t(end)/60]);

% density and mass in storage, figure 4
figure;
subplot(2,2,1);
plot(data.t/60,data.rhov1,data.t/60,data.rhov2,data.t/60,data.ETTTVenstate)
ylabel('Vapor densities (g/L)');
legend('(ST)','(ET)','(ET) vent state');
xlim([0 data.t(end)/60]);
xlabel('Time (min)');
grid on;

subplot(2,2,2);
plot(data.t/60,data.rho_L1,data.t/60,data.rho_L2)
ylabel('Liquid densities (g/L)');
legend('(ST)','(ET)');
xlim([0 data.t(end)/60]);
xlabel('Time (min)');
grid on;

subplot(2,2,3:4);
plot(data.t/60,data.mL2+data.mv2,data.t/60,data.mL2);
ylabel('Total mass in ET (kg)');
xlim([0 data.t(end)/60]);
legend('Total mass','liquid only');
xlabel('Time (min)');
grid on;

% venting fluxes, figure 5
figure;
subplot(2,2,1:2);
plot(data.t/60,data.Boiloff_ET)
ylabel('Transfer losses from (ET) (kg)');
xlabel('Time (min)');
legend('Transfer losses');
grid on;
xlim([0 data.t(end)/60]);

subplot(2,2,3:4);
plot(data.t/60,data.mL1+data.mv1,data.t/60,data.mL1)
ylabel('Mass in (ST) (kg)');
xlim([0 data.t(end)/60]);
xlabel('Time (min)');
legend('total mass','liquid only');
grid on;

% mass flows for vapor phases, figure 6
figure
subplot(2,2,1:2);
%AAA = Jvvalve1;
plot(data.t/60,data.Jv10,data.t/60,data.Jcd1,data.t/60,data.AAA,data.t/60,data.Jboil);
ylabel('mass flows in VAPOR (ST), in kg/sec');
xlim([0 data.t(end)/60]);
xlabel('Time (min)');
legend('Jv1', 'Jcd1','Jvvalve1','Jboil');
grid on;

subplot(2,2,3:4);
plot(data.t/60,data.Jv20,data.t/60,data.Jcd2,data.t/60,data.Jvalve222);
ylabel('mass flows in VAPOR (ET), in kg/sec');
xlim([0 data.t(end)/60]);
xlabel('Time (min)');
legend('Jv2', 'Jcd2','Jvvalve2');
grid on;

% energy balance in (ST) and (ET), figure 7
figure;
subplot(2,2,1);
plot(data.t/60,(40)*ones(size(data.t)),data.t/60,data.BBB,data.t/60,-data.CCC,data.t/60,data.DDD,data.t/60,data.EEE,data.t/60,data.KKK,data.t/60,data.FFF)
legend('QdotEV',' - QdotVS','-pdV','-Jvalve.h','-Jcd.h','QdotV','Jboil.h' )
ylabel('Heat transfers in VAPOR - ST (Watts)');
xlim([0 data.t(end)/60]);
xlabel('Time (min)');
grid on;

subplot(2,2,2);
plot(data.t/60,(200)*ones(size(data.t)),data.t/60,data.GGG,data.t/60,data.CCC,data.t/60,data.HHH,data.t/60,data.III,data.t/60,data.LLL,data.t/60,data.JJJ)
legend('QdotEL',' - QdotLS','pdV','-Jtr.h','Jcd.h','QdotL','-Jvap.h' )
xlim([0 data.t(end)/60]);
xlabel('Time (min)');
ylabel('Heat transfers LIQUID - ST (Watts)');
grid on;

subplot(2,2,3);
plot(data.t/60,data.MMM,data.t/60,data.NNN,data.t/60,-data.OOO,data.t/60,data.QQQ,data.t/60,data.RRR,data.t/60,data.WWW,data.t/60,data.PPP)
legend('QdotWV',' - QdotVS','-pdV','-Jvalve.h','-Jcd.h','QdotV','topfill' )
ylabel('Heat transfers in VAPOR - ET (Watts)');
xlim([0 data.t(end)/60]);
xlabel('Time (min)');
grid on;

subplot(2,2,4);
plot(data.t/60,data.SSS,data.t/60,data.TTT,data.t/60,data.OOO,data.t/60,data.UUU,data.t/60,data.VVV,data.t/60,data.XXX)
legend('QdotWL',' - QdotLS','pdV','Jtr.h','Jcd.h','QdotL' )
xlim([0 data.t(end)/60]);
xlabel('Time (min)');
ylabel('Heat transfers LIQUID - ET (Watts)');
grid on;



