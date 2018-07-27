% Script use to run the LH2 simulation code

% here : (ST) or 1 is trailer (horizontal cylinder)
%        (ET) or 2 is station storage (vertical cylinder)   
clc
close all
tic;
% 1. initialize parameters
inputs_TrailerToDewar;

% 2. run simulation with default timespan
nominal= LH2Simulate;

% 3. extract data
display('now extracting data');
Data_extraction;

% 4. plot results
plotLH2Data(nominal);

% 5. save data
dlmwrite('output_bottomfill.txt',cat(2,nominal.t(1:10:end),nominal.t(1:10:end)/60, ...
    nominal.mL2(1:10:end), nominal.mv2(1:10:end),nominal.mL2(1:10:end)+nominal.mv2(1:10:end),nominal.Boiloff_ET(1:10:end),...
    nominal.TL2((1:10:end),:),nominal.Ts2((1:10:end),:),nominal.Tv2((1:10:end),:),nominal.Tw2((1:10:end),:),...
    nominal.uL2((1:10:end),:),nominal.uv2((1:10:end),:), ...
    nominal.rho_L2(1:10:end), (nominal.rhov2(1:10:end)), ...    
    nominal.pL2(1:10:end)/6894.75729, nominal.pv2(1:10:end)/6894.75729,(nominal.hL2(1:10:end)), ...
    nominal.Jtr(1:10:end),nominal.Jvalve222(1:10:end),nominal.Jcd2(1:10:end),...
    nominal.MMM(1:10:end),nominal.NNN(1:10:end),nominal.OOO(1:10:end),nominal.PPP(1:10:end),nominal.QQQ(1:10:end),...
    nominal.RRR(1:10:end),nominal.SSS(1:10:end),nominal.TTT(1:10:end),...
    nominal.UUU(1:10:end),nominal.VVV(1:10:end),nominal.WWW(1:10:end),nominal.XXX(1:10:end),...
    nominal.mL1(1:10:end), nominal.mv1(1:10:end),nominal.mL1(1:10:end)+ nominal.mv1(1:10:end),...
    nominal.TL1((1:10:end),:),nominal.Ts1((1:10:end),:),nominal.Tv1((1:10:end),:), ...
    nominal.uL1((1:10:end),:),nominal.uv1((1:10:end),:), ...
    nominal.rho_L1(1:10:end), (nominal.rhov1(1:10:end)), ...
    nominal.pL1(1:10:end)/6894.75729, nominal.pv1(1:10:end)/6894.75729,(nominal.hL1(1:10:end)'), ...
    nominal.AAA(1:10:end), nominal.Jcd1(1:10:end), nominal.Jboil(1:10:end),...
    nominal.BBB(1:10:end),nominal.CCC(1:10:end),nominal.DDD(1:10:end),nominal.EEE(1:10:end),nominal.FFF(1:10:end),...
    nominal.GGG(1:10:end),nominal.HHH(1:10:end),nominal.III(1:10:end),...
    nominal.JJJ(1:10:end),nominal.KKK(1:10:end),nominal.LLL(1:10:end)...
    ));
toc;
    %AAA = Jvvalve1;
    %BBB = - QdotVS1;
    %CCC = pdV1;
    %DDD = - Jvvalve1 *(hvalve1+0.5*vv1^2-((refpropm('U','T',Tv1(P.nV1),'D',rhov1,'PARAHYD')+P.refstatedelta)));
    %EEE =  - Jcd1*(hcd1-((refpropm('U','T',Tv1(P.nV1),'D',rhov1,'PARAHYD')+P.refstatedelta)));
    %FFF = Jboil*(hboil-((refpropm('U','T',Tv1(P.nV1),'D',rhov1,'PARAHYD')+P.refstatedelta)));
    %GGG = - QdotLS1;
    %HHH = - Jtr*(htr+0.5*vtr^2- (refpropm('U','T',TL1(P.nL1),'Q',0,'PARAHYD')+P.refstatedelta));
    %III = + Jcd1*(hcd - (refpropm('U','T',TL1(P.nL1),'Q',0,'PARAHYD')+P.refstatedelta));
    %JJJ = - Jvap*(htr - (refpropm('U','T',TL1(P.nL1),'Q',0,'PARAHYD')+P.refstatedelta));
    %KKK = QdotV1 - Jv1*(refpropm('U','T',Tv1(P.nV1),'D',rhov1,'PARAHYD')+P.refstatedelta);
    %LLL = QdotL1 - JL1*(refpropm('U','T',TL1(P.nL1),'Q',0,'PARAHYD')+P.refstatedelta);
    
    %MMM = QdotWV2;
    %NNN = - QdotVS2;
    %OOO = pdV2;   
    %PPP = P.ratio_top_bottom*Jtr*(htr+0.5*vtr^2-(refpropm('U','T',Tv2(P.nV2),'D',rhov2,'PARAHYD')+P.refstatedelta));
    %QQQ = - Jvvalve2*(hvalve2 + 0.5*vv2^2-(refpropm('U','T',Tv2(P.nV2),'D',rhov2,'PARAHYD')+P.refstatedelta));
    %RRR = - Jcd2*(hcd2-(refpropm('U','T',Tv2(P.nV2),'D',rhov2,'PARAHYD')+P.refstatedelta));
    %SSS = QdotWL2;
    %TTT = - QdotLS2;
    %UUU = + (1-P.ratio_top_bottom)*Jtr*(htr+0.5*vtr^2-(refpropm('U','T',TL2(P.nL2),'D',rho_L2,'PARAHYD')+P.refstatedelta));
    %VVV =  + Jcd2*(hcd2-(refpropm('U','T',TL2(P.nL2),'D',rho_L2,'PARAHYD')+P.refstatedelta));
    %WWW = QdotV2 - Jv2*(refpropm('U','T',Tv2(P.nV2),'D',rhov2,'PARAHYD')+P.refstatedelta);
    %XXX = QdotL2 - JL2*(refpropm('U','T',TL2(P.nL2),'D',rho_L2,'PARAHYD')+P.refstatedelta);  