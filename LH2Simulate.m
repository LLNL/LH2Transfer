function data = LH2Simulate(name)
% main code for LH2 transfer simulation

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

% set default name
if nargin<2
	name = 'fill from trailer to Dewar';
end

P.waitbar = waitbar(0,['Simulating ' name '...']);

% set up initial state
UL10 = refpropm('U','T',P.TL10,'Q',0,'PARAHYD')*ones(P.nL1,1);
Uv10 = refpropm('U','T',P.Tv10,'Q',1,'PARAHYD')*ones(P.nV1,1);
UL20 = refpropm('U','T',P.TL20,'Q',0,'PARAHYD')*ones(P.nL2,1);
Uv20 = refpropm('U','T',P.Tv20,'Q',1,'PARAHYD')*ones(P.nV2,1);

x0 = [  P.mL10;
		UL10;
		P.mv10;
        Uv10;
		P.Ts10;
		P.Jtr0;
		P.mVap0;
		P.Jboil0;
		P.mL20;
        UL20;
		P.mv20;
        Uv20;
        P.Ts20;
		P.Tw20;
        0;0;0;0;
        0;0;0;0;
        0;0;0;0;
        0;0;0;0;
        0;0;0;0;
        0;0;0;0;
        0;0
];

% declare and initialize global variables
global ETTVentState;    % state of venting valve for (ET), 0 or 1
ETTVentState = P.ETVentState; % initial value
global ET_fill_complete; % flag indicating (ET) is completely full, per LH2Model.TopET criteria
ET_fill_complete = 0;    % initial value for flag
global ST_ready;         % flag inidicating (ST) is ready to deliver fuel, i.e. pv1 has reached delivery pressure
ST_ready = 0;            % initial value for flag
global ST_vent_complete; % flag indicating (ST) vent is complete, per LH2Model.p_ST_final criteria
ST_vent_complete = 0;    % initial value for flag

% initialize variables for ODE solver
tstart = 0;
tout = tstart;
xout = x0';
teout = [];
xeout = [];
ieout = [];
tfinal=P.tFinal;
xout=horzcat(xout,ETTVentState);

function dxdt = LH2dxdt(P,t,x)
                
    display('time, in min');
    display(t/60);

    %waitbar
    waitbar(t/P.tFinal,P.waitbar);

    % obtain initial state variables
    mL1 = x(1);
    uL1 = x(2:1+P.nL1);
    mv1 = x(P.nL1+2);
    uv1 = x(P.nL1+3:P.nL1+2+P.nV1);
    Ts1 = x(P.nL1+P.nV1+3);
    Jtr = x(P.nL1+P.nV1+4);
    mVap = x(P.nL1+P.nV1+5);
    Jboil = x(P.nL1+P.nV1+6);
    mL2 = x(P.nL1+P.nV1+7);
    uL2 = x(P.nL1+P.nV1+8:P.nL1+P.nV1+P.nL2+7);
    mv2 = x(P.nL1+P.nV1+P.nL2+8);
    uv2 = x(P.nL1+P.nV1+P.nL2+9:P.nL1+P.nV1+P.nL2+P.nV2+8);
    Ts2 = x(P.nL1+P.nV1+P.nL2+P.nV2+9);
    Tw2 = x(P.nL1+P.nV1+P.nL2+P.nV2+10);
    
    % to make sure Ts is not too low (<14 K have been observed under very stiff conditions)
    if Ts1<14 
        Ts1=14;
    end
    
    if Ts2<14
        Ts2=14;
    end
    
    % raises error if mL1 <= 0
    if mL1 <=0
        error('liquid mass equal to zero in (ST)');
    end
      
    %---------------------
    % ST initial calculations
    %---------------------
    
    rho_L1 = -5.12074746E-07*(uL1(P.nL1)/1000)^3 - 1.56628367E-05*(uL1(P.nL1)/1000)^2 - 1.18436797E-01*(uL1(P.nL1)/1000) + 7.06218354E+01; % correlation from REFPROP v9.1
    VL1 = mL1/rho_L1;               % [m^3] volume of liquid in (ST)
    hL1 = cylVToH(VL1,P.R1,P.Lcyl); % [m] height of liquid in (ST)              
    Vullage1 = P.VTotal1-VL1;       % [m^3] volume of vapor in (ST)
    rhov1 = mv1/Vullage1;           % [g/l] density of vapor in (ST)
    
    try
      Tv1(P.nV1) = refpropm('T','D',rhov1,'U',uv1(P.nV1),'PARAHYD'); % temperature of the bulk vapor in (ST)
    catch
      Tv1(P.nV1) = refpropm('T','D',rhov1,'U',fix(100*uv1(P.nV1))/100,'PARAHYD'); % temperature of the bulk vapor in (ST)
    end

    Tv1(P.nV1)= fix(Tv1(P.nV1)*100)/100; % rounding to the third decimal. This is to avoid non-convergence error around the critical point for REPRPOP. rounding() does not seem to work adequately.

    pv1=vaporpressure(uv1(P.nV1),rhov1); % [Pa] vapor pressure in (ST)
    pL1 = rho_L1*P.g*hL1;                % [Pa] pressure at bottom of (ST) dues to liquid weight
    pTotal1 = pv1+pL1;                   % [Pa] total pressure in (ST)
    
    % this is to make sure that the pressure in (ST) is large enough for delivery
    if pv1 > min(P.p_ST_slow,P.p_ST_fast)
        ST_ready = 1;
    end
    
    % this is to check whether vent from (ST) is complete
    if pv1 <= P.p_ST_final && ET_fill_complete == 1
        ST_vent_complete = 1;
    end
    
    % liquid temperatures for (ST)
    for i = 1:P.nL1
        TL1(i)= 1.44867559E-07*(uL1(i)/1000)^3 - 2.53438808E-04*(uL1(i)/1000)^2 + 1.05449468E-01*(uL1(i)/1000) + 2.03423757E+01; % correlation from REFPROP v9.1
        if TL1(i) < 13.804 
            TL1(i) = 13.804;
        elseif TL1(i) > 32.93
            TL1(i) = 32.93;
        end
   end
    
    % vapor temperatures for (ST)
    for i = 1:P.nV1-1
        if uv1(i)<0
            uv1(i)=0;
        end
        display(i);
        display(pv1/1000);
        display(uv1(i));
       Tv1(i) = refpropm('T','P',pv1/1000,'U',uv1(i)/1.5,'PARAHYD'); % temperature of the vapor layers  in (ST)
    end
       
      
    % computation of the surface area between vapor and liquid in (ST), that is a horizontal cylinder
    if hL1 > P.R1
        d = hL1 - P.R1;
    else
        d = P.R1 - hL1;
    end
    c = 2 * P.R1 * sqrt(1-(d/P.R1)^2);
    S1 = c * P.Lcyl; % new area of interface between vapor and liquid
    
    
    %----------------------    
    % ET initial calculations
    %---------------------
    rho_L2 = -5.12074746E-07*(uL2(P.nL2)/1000)^3 - 1.56628367E-05*(uL2(P.nL2)/1000)^2 - 1.18436797E-01*(uL2(P.nL2)/1000) + 7.06218354E+01; % correlation from REFPROP v9.1
    VL2 = mL2/rho_L2;         % [m^3] volume of liquid in (ET)
    Vullage2 = P.VTotal2-VL2; % [m^3] ullage volume in (ET)
    rhov2 = mv2/Vullage2;     % [g/L] vapor density in (ET)   
   
   try
     quality2=refpropm('q','D',rhov2,'U',uv2(P.nV2),'PARAHYD');
   catch
     display('Non-truncation of uv2 did not converge in "quality2". Choosing truncated value instead');
     uv2(P.nV2) = fix(100*uv2(P.nV2))/100;
     %display(uv2(P.nV2));
     try
        quality2=refpropm('q','D',rhov2,'U',uv2(P.nV2),'PARAHYD');
     catch
         quality2=refpropm('q','D',fix(100*rhov2)/100,'U',uv2(P.nV2),'PARAHYD');
         %uv2(P.nV2) = fix(100*uv2(P.nV2))/100;
         rhov2 = fix(100*rhov2)/100;
     end
   end
   
   try
       pv2=vaporpressure(uv2(P.nV2),rhov2); % [Pa] vapor presusre in (ET)
   catch
       display('Non-truncation of uv2 did not converge in "Pv2". Choosing truncated value instead');
       pv2=vaporpressure(fix(100*uv2(P.nV2))/100,rhov2); % [Pa] vapor presusre in (ET)
   end
   
   if pv2==0
       display('!!! Pv2=0 !!!'); 
       pv2=vaporpressure(fix(100*uv2(P.nV2))/100,rhov2);
       uv2(P.nV2)= fix(100*uv2(P.nV2))/100;        
   end
   
   try
        Tv2(P.nV2) = refpropm('T','D',rhov2,'U',uv2(P.nV2),'PARAHYD'); % [K] bulk vapor temperature in (ET)\ display('truncation of uv2 did not converge in "quality2". Choosing non-truncated value instead');
    catch
        Tv2(P.nV2) = refpropm('T','D',rhov2,'U',fix(100*uv2(P.nV2))/100,'PARAHYD'); % [K] bulk vapor temperature in (ET)
        display('Non-truncation of uv2 did not converge in "Tv2". Choosing truncated value instead');
    end
    
    Tv2(P.nV2)= fix(Tv2(P.nV2)*100)/100; % rounding to the third decimal. This is to avoid non-convergence error around the critical point for REPRPOP. rounding() does not seem to work adequately.

    if quality2 > 0 && quality2 < 1 %
        liquiddensity=-5.24588E-05*Tv2(P.nV2)^6 + 7.39502E-03*Tv2(P.nV2)^5 - 4.29976E-01*Tv2(P.nV2)^4 + 1.31922E+01*Tv2(P.nV2)^3 - 2.25208E+02*Tv2(P.nV2)^2 + 2.02705E+03*Tv2(P.nV2) - 7.43508E+03;
        hL2 = (VL2+(1-quality2)*mv2/liquiddensity)/P.A2;         % [m] height of liquid in (ET)
    else
        hL2 = VL2/P.A2;         % [m] height of liquid in (ET)
    end
     
    if hL2 > P.TopET * P.H  % stopping criteria for (ET) filling
        ET_fill_complete = 1;
    end

    pL2 = rho_L2*P.g*hL2; % [Pa] liquid pressure in (ET)
    pTotal2 = pv2+pL2;    % [Pa] total pressure in (ET)
  
    % liquid temperature in (ET)
    for i = 1:P.nL2
         TL2(i)= 1.44867559E-07*(uL2(i)/1000)^3 - 2.53438808E-04*(uL2(i)/1000)^2 + 1.05449468E-01*(uL2(i)/1000) + 2.03423757E+01; % correlation from REFPROP v9.1
        if TL2(i) < 13.804 
            display('Low TL2');
            display(TL2(i));
            TL2(i) = 13.804;
        elseif TL2(i) > 32.93
            display('High TL2');
            display(TL2(i));
            TL2(i) = 32.93;
        end 
    end
    
    % vapor temperatures in (ET)
    for i = 1:P.nV2-1
        if uv2(i)<0
            uv2(i)=0;
        end
        display(i);
        display(pv2/1000);
        display(uv2(i));
      % Tv2(i) = refpropm('T','U',uv2(i),'P',pv2/1000,'PARAHYD');
       Tv2(i) = refpropm('T','P',pv2/1000,'U',uv2(i),'PARAHYD');
       display(Tv2(i));
    end

    %---------------------
    % surface temperatures
    %---------------------
    Ts10 = P.T_c*(pv1/P.p_c)^(1/P.lambda); % From Osipov 2008, see reference in Readme file
    Ts20 = P.T_c*(pv2/P.p_c)^(1/P.lambda);   

    dTs1dt = (Ts10-Ts1)/P.tminL1;
    dTs2dt = (Ts20-Ts2)/P.tminL2;

    % enthalpy of vaporization % correlation from REFPROP v9.1
    qh1 = 1000 * (-0.002445451720487*Ts1^6 + 0.3629946692976*Ts1^5 - 22.28028769483*Ts1^4 + 723.6541112107*Ts1^3 - 13116.31006512*Ts1^2 + 125780.2915522*Ts1- 498095.5392318);
    qh2 = 1000 * (-0.002445451720487*Ts2^6 + 0.3629946692976*Ts2^5 - 22.28028769483*Ts2^4 + 723.6541112107*Ts2^3 - 13116.31006512*Ts2^2 + 125780.2915522*Ts2- 498095.5392318);

    % obtain control inputs
    U = LH2Control(hL2,pv1,pv2, ET_fill_complete,ST_vent_complete);

    % calculate transmission line parameters
    apipe = 2*pi*(P.DPipe/2)^2*sqrt(rho_L1*P.DPipe/2/P.LPipe/P.f);
    
    % calculate valve area and lambda for fill valve
    AE = (2*pi*(P.dE/2)^2);
    lambdaE = U.lambdaE;
    alphaE = AE*sqrt(2*rho_L1/P.kE);
    if lambdaE <= 0 % this is to avoid divided by 0 error
        aeff = 0;
    else
         aeff =  ((lambdaE*alphaE)^-2 + apipe^-2 )^(-1/2);
    end
    
    Jtr0 = ST_ready*aeff*dsqrt(pTotal1-pTotal2);
    dJtrdt = (Jtr0 - Jtr)/P.tau_tr;
      
    %---------------------
    % vaporizer in (ST)
    %---------------------
    Jvap = abs(P.c_vap*U.lambdaV*dsqrt(2*rho_L1*(pTotal1-P.p_atm))); % flow into vaporizer. abs() added so that small quantities in ST can work.

    if mVap<=0
        Jboil0 = 0;
        Jboil = max(0,Jboil);
        dmVapdt = max(0,Jvap - Jboil);
    else
        Jboil0 = Jvap;
        dmVapdt = Jvap - Jboil; % vaporizer mass flow
    end
    dJboildt = (Jboil0 - Jboil)/P.tau_vap;
    
    P.VapValveState = U.lambdaV;

    % determine ST vent valve state
    P.STVentState = U.STVentState;% store STVentState value for next iteration
    
    %----------------------
    % vent flows
    %----------------------
    % compute vapor flow for ST end ET vents
    Jvvalve1 = P.STVentState*gasFlow(P.S_valve1,P.gamma_,rhov1,pv1,P.p_atm); 
    Jvvalve2 = ETTVentState * gasFlow(P.S_valve2,P.gamma_,rhov2,pv2,P.p_atm);
    
    %---------------------
    % heat transfer between the saturated film and the vapor and liquid phases in (ST)
    %---------------------
    % transport properties in the gas and liquid phases
    if refpropm('Q','D',rhov1,'U',uv1(P.nV1),'PARAHYD') < 1 && refpropm('Q','D',rhov1,'U',uv1(P.nV1),'PARAHYD') > 0     
        P.kappa_v = refpropm('L','T',Tv1(P.nV1),'Q',1,'PARAHYD'); % REFPROP does not work here for qualities different than 1 or 0
        P.mu_v = refpropm('V','T',Tv1(P.nV1),'Q',1,'PARAHYD');
        P.cv_v = refpropm('O','T',Tv1(P.nV1),'Q',1,'PARAHYD'); % Cv
        P.cp_v = refpropm('C','T',Tv1(P.nV1),'Q',1,'PARAHYD'); % Cp
        beta_v = refpropm('B','T',Tv1(P.nV1),'Q',1,'PARAHYD');
    else
        P.kappa_v = refpropm('L','D',rhov1,'U',uv1(P.nV1),'PARAHYD');
        P.mu_v = refpropm('V','D',rhov1,'U',uv1(P.nV1),'PARAHYD');
        P.cv_v = refpropm('O','D',rhov1,'U',uv1(P.nV1),'PARAHYD');% Cv
        P.cp_v = refpropm('C','D',rhov1,'U',uv1(P.nV1),'PARAHYD');% Cp
        beta_v = refpropm('B','D',rhov1,'U',uv1(P.nV1),'PARAHYD');
    end
    
     P.kappa_L = refpropm('L','T',TL2(P.nL1),'Q',0,'PARAHYD');
     P.cv_L = refpropm('O','T',TL1(P.nL1),'Q',0,'PARAHYD');% Cv
     P.cp_L = refpropm('C','T',TL1(P.nL1),'Q',0,'PARAHYD');% Cp
     beta_L = refpropm('B','T',TL1(P.nL1),'Q',0,'PARAHYD');
     P.mu_L = refpropm('V','T',TL1(P.nL1),'Q',0,'PARAHYD');
      
    % set up grid for vapor in ST (count from interface to top)
    lmin = sqrt(P.kappa_v*P.tminV1/P.c_v/rhov1);
    l_V1(1) = lmin/(1+exp(pi/2/sqrt(P.nV1)));			% h_0
    l12_V1(1) = lmin;									% h_1/2
    for i=2:P.nV1
        l12_V1(i) = l12_V1(i-1)*exp(pi/sqrt(P.nV1));	% h_i+1/2
        l_V1(i) = sqrt(l12_V1(i-1)*l12_V1(i));			% h_i
    end
    
    % set up grid for liquid in ST (count from interface to bottom)
    lmin = sqrt(P.kappa_L*P.tminL1/P.cv_L/rho_L1);
    l_L1(1) = lmin/(1+exp(pi/2/sqrt(P.nL1)));			% h_0
    l12_L1(1) = lmin;									% h_1/2
    for i=2:P.nL1
        l12_L1(i) = l12_L1(i-1)*exp(pi/sqrt(P.nL1));	% h_i+1/2
        l_L1(i) = sqrt(l12_L1(i-1)*l12_L1(i));			% h_i
    end  
   
    hVS1_cond = P.kappa_v/l12_V1(1);
    hVS1_conv = P.kappa_v*0.156*(P.g*beta_v*P.cp_v*rhov1^2*(Ts1-Tv1(P.nV1))/P.kappa_v/P.mu_v)^(1/3);
    
    hLS1_cond = P.kappa_L/l12_L1(1);
    hLS1_conv = P.kappa_L*0.156*(P.g*beta_L*P.cp_L*rho_L1^2*abs(TL1(P.nL1)-Ts1)/P.kappa_L/P.mu_L)^(1/3);
        
    % heat flow terms (ST)
    QdotLS1_cond = hLS1_cond*S1*(TL1(1)-Ts1) - l_L1(1)*P.cp_L*rho_L1*dTs1dt;
    QdotLS1_conv = hLS1_conv*S1*(TL1(1)-Ts1)*(TL1(P.nL2)>Ts1);
    if QdotLS1_conv>0
        QdotLS1 = max(QdotLS1_conv,QdotLS1_cond);
    else
        QdotLS1 = QdotLS1_cond; % Q_dotLS1_conv is 0 here
    end

    QdotVS1_conv = hVS1_conv*S1*(Tv1(1)-Ts1)*(Ts1>Tv1(P.nV1));
    QdotVS1_cond = hVS1_cond*S1*(Tv1(1)-Ts1) - l_V1(1)*P.cv_v*rhov1*dTs1dt;
    if QdotVS1_conv<0
        QdotVS1 = min(QdotVS1_conv,QdotVS1_cond);
    else
        QdotVS1 = QdotVS1_cond; % Q_dotVS1_conv is 0 here
    end

    %---------------------
    % heat transfer between the wall and the vapor and liquid phases in (ET)
    %---------------------
    % transport properties for convection between wall and vapor
    if quality2 < 1 && quality2 > 0
    %if refpropm('Q','D',rhov2,'U',uv2(P.nV2),'PARAHYD') < 1 && refpropm('Q','D',rhov2,'U',uv2(P.nV2),'PARAHYD') > 0
        Pr = refpropm('^','T',Tv2(P.nV2),'Q',1,'PARAHYD'); % REFPROP does not work here for qualities different than 1 or 0
        P.kappa_v = refpropm('L','T',Tv2(P.nV2),'Q',1,'PARAHYD');
        P.mu_v = refpropm('V','T',Tv2(P.nV2),'Q',1,'PARAHYD');
        P.cv_v = refpropm('O','T',Tv2(P.nV2),'Q',1,'PARAHYD'); % Cv
        P.cp_v = refpropm('C','T',Tv2(P.nV2),'Q',1,'PARAHYD'); % Cp
        beta_v = refpropm('B','T',Tv2(P.nV2),'Q',1,'PARAHYD');
    else
        Pr = refpropm('^','D',rhov2,'U',uv2(P.nV2),'PARAHYD');
        P.kappa_v = refpropm('L','D',rhov2,'U',uv2(P.nV2),'PARAHYD');
        P.mu_v = refpropm('V','D',rhov2,'U',uv2(P.nV2),'PARAHYD');
        P.cv_v = refpropm('O','D',rhov2,'U',uv2(P.nV2),'PARAHYD');% Cv
        P.cp_v = refpropm('C','D',rhov2,'U',uv2(P.nV2),'PARAHYD');% Cp
        beta_v = refpropm('B','D',rhov2,'U',uv2(P.nV2),'PARAHYD');
    end
    
    nuv2 = P.mu_v/rhov2;  % viscosity for vapor phase
    Ra = abs(P.g* beta_v*(Tw2-Tv2(P.nV2))*(P.H-hL2)^3*Pr/nuv2);
    Psi = (1+(0.492/Pr)^(9/16))^(-16/9);
    Nu = 0.68+0.503*(Ra*Psi)^(1/4);

    % transport properties for convection between wall and liquid
    Pr_L = refpropm('^','T',TL2(P.nL2),'Q',0,'PARAHYD');
    P.kappa_L = refpropm('L','T',TL2(P.nL2),'Q',0,'PARAHYD');
    P.mu_L = refpropm('V','T',TL2(P.nL2),'Q',0,'PARAHYD');
    P.cv_L = refpropm('O','T',TL2(P.nL2),'Q',0,'PARAHYD');% Cv
    P.cp_L = refpropm('C','T',TL2(P.nL2),'Q',0,'PARAHYD');% Cp
    beta_L = refpropm('B','T',TL2(P.nL2),'Q',0,'PARAHYD');
    nuL2 = P.mu_L/rho_L2; % viscosity for liquid phase
      
    Ra_L_side = abs(P.g*beta_L*(Tw2-TL2(P.nL2))*(hL2)^3*Pr_L/nuL2);
    Psi_L_side = (1+(0.492/Pr_L)^(9/16))^(-16/9); % Nusselt correlation for flow along vertical wall
    Nu_L_side = 0.68+0.503*(Ra_L_side*Psi_L_side)^(1/4);% Nusselt correlation for flow along vertical wall
    Ra_L_bottom = abs(P.g*beta_L*(Tw2-TL2(P.nL2))*(P.R2/2)^3*Pr_L/nuL2); 
    Nu_L_bottom = 0.27 * (Ra_L_bottom)^(1/4);% Nusselt correlation for flow along horizontal plate
    
    
    hWV2 = Nu*P.kappa_v/(P.H-hL2);
    hWL2_side = Nu_L_side*P.kappa_L/(hL2);  
    hWL2_bottom = Nu_L_bottom*P.kappa_L/(P.R2/2);

    QdotWL2 = (Tw2-TL2(P.nL2))*(hWL2_bottom*P.A2+hWL2_side*2*pi*P.R2*hL2); % convection term, wall to liquid
    QdotWV2 = hWV2*(Tw2-Tv2(P.nV2))*(P.A2+2*pi*P.R2*(P.H-hL2));            % convection term, wall to vapor
    
    %---------------------       
    % heat transfer between the saturated film and the vapor / liquid phases in (ET)
    %---------------------
    
    % set up grid for liquid in ET (count from interface to bottom)
    lmin = sqrt(P.kappa_L*P.tminL2/P.cv_L/rho_L2);
    l_L2(1) = lmin/(1+exp(pi/2/sqrt(P.nL2)));			% h_0
    l12_L2(1) = lmin;									% h_1/2
    for i=2:P.nL2
        l12_L2(i) = l12_L2(i-1)*exp(pi/sqrt(P.nL2));	% h_i+1/2
        l_L2(i) = sqrt(l12_L2(i-1)*l12_L2(i));			% h_i
    end
    
    % set up grid for vapor in ET  (count from interface to top)
    lmin = sqrt(P.kappa_v*P.tminV2/P.cv_v/rhov2);
    l_V2(1) = lmin/(1+exp(pi/2/sqrt(P.nV2)));			% h_0
    l12_V2(1) = lmin;									% h_1/2
    for i=2:P.nV2
        l12_V2(i) = l12_V2(i-1)*exp(pi/sqrt(P.nV2));	% h_i+1/2
        l_V2(i) = sqrt(l12_V2(i-1)*l12_V2(i));			% h_i
    end
    
    hVS2_cond = P.kappa_v/l12_V2(1);
    hVS2_conv = P.kappa_v*0.156*(P.g*beta_v*P.cp_v*rhov2^2*(Ts2-Tv2(P.nV2))/P.kappa_v/P.mu_v)^(1/3);
    hLS2_cond = P.kappa_L/l12_L2(1);
    hLS2_conv = P.kappa_L*0.156*(P.g*beta_L*P.cp_L*rho_L2^2*abs(TL2(P.nL2)-Ts2)/P.kappa_L/P.mu_L)^(1/3);

    QdotLS2_cond = hLS2_cond*P.A2*(TL2(1)-Ts2) - l_L2(1)*P.cp_L*rho_L2*dTs2dt;
    QdotLS2_conv = hLS2_conv*P.A2*(TL2(1)-Ts2)*(TL2(P.nL2)>Ts2);
    
    if QdotLS2_conv>0
        QdotLS2 = max(QdotLS2_conv,QdotLS2_cond);
    else
        QdotLS2 = QdotLS2_cond; % Q_dotLS2_conv is 0 here
    end

    QdotVS2_conv = hVS2_conv*P.A2*(Tv2(1)-Ts2)*(Ts2>Tv2(P.nV2));
    QdotVS2_cond = hVS2_cond*P.A2*(Tv2(1)-Ts2) - l_V2(1)*P.cv_v*rhov2*dTs2dt;
    
    if QdotVS2_conv<0
        QdotVS2 = min(QdotVS2_conv,QdotVS2_cond);
    else
        QdotVS2 = QdotVS2_cond; % Q_dotVS2_conv is 0 here
    end

    %---------------------    
    % condensation flows (ST) and (ET)
    %---------------------
    if qh1 < 0
        Jcd1 = 0;
    else
        Jcd1 = -(QdotLS1+QdotVS1)/qh1;
    end
    
    if qh2<=0
        Jcd2 = 0;
    else       
        Jcd2 = -(QdotLS2+QdotVS2)/qh2 - 0*(P.ratio_top_bottom) * Jtr; % term added for top fill
    end
      
    %---------------------    
    % mass balances (ST) and (ET)
    %---------------------
    Jv1 = Jboil - Jvvalve1 - Jcd1;                      % variation of mass of vapor in (ST)
    JL1 = -Jtr - Jvap + Jcd1;                           % variation of mass of liquid in (ST)
        
    Jv2 = (P.ratio_top_bottom) * Jtr  - Jvvalve2 - Jcd2; % variation of mass of vapor in (ET)
    JL2 = (1-P.ratio_top_bottom) * Jtr + Jcd2;         % variation of mass of liquid in (ET)

    %---------------------
    % pdV work, ST and ET
    %---------------------
    pdV1 = -pv1*(JL1/rho_L1);
    pdV2 = -pv2*(JL2/rho_L2);
    
    %---------------------
    % exit velocities, ST and ET
    %---------------------
    vv1 = Jvvalve1/P.S_valve1/rhov1;
    vv2 = Jvvalve2/P.S_valve2/rhov2;
    
    %---------------------
    % enthalpy terms, modified for ideal vs. real gases
    %---------------------
    if TL1(P.nL1) > 32 
        htr_L = P.c_L*TL1(P.nL1); 
    else
        if TL1(P.nL1) < 14
            TL1(P.nL1) = 14;
        end  
        htr_L= refpropm('H','T',TL1(P.nL1),'Q',0,'PARAHYD');
    end
    
    if Ts1 > 32
         hcd1 = P.c_p*Ts1;
    else
         hcd1 = refpropm('H','T',Ts1,'Q',1,'PARAHYD');
    end   

    if Ts2 > 32
         hcd2 = P.c_p*Ts2;
    else
        hcd2 = refpropm('H','T',Ts2,'Q',1,'PARAHYD');
    end
    
    if P.Tboil > 32 % Tboil is the assumed temperature of the boiling molecules
        hboil = P.c_p*P.Tboil;
    else 
        hboil = refpropm('H','T',P.Tboil,'Q',1,'PARAHYD');
    end
    
    hvalve1 = refpropm('H','T',Tv1(P.nV1),'D',rhov1,'PARAHYD');
    hvalve2 = refpropm('H','T',Tv2(P.nV2),'D',rhov2,'PARAHYD');
      
    %------------------------------------------------------
    % Heat flows to vapor and liquid phases in (ST)
    %------------------------------------------------------
    % heat flow to vapor in (ST)
    QdotV1 = P.QdotEV1 - QdotVS1 - pdV1 ...
        - Jvvalve1 *(hvalve1+0.5*vv1^2) ...
        - Jcd1*hcd1 ...
        + Jboil*hboil;
        
    % heat flow to liquid in (ST)
    rhotr = rho_L1 ;                    % assumed density in the transfer line
    vtr = Jtr/(pi*(0.5*P.dE)^2)/rhotr;  % velocity in the transfer line
    
    QdotL1 = P.QdotEL1 - QdotLS1 + pdV1 ... 
            - Jtr*(htr_L+0.5*vtr^2) ...
            + Jcd1*hcd1...
            - Jvap*htr_L;

   %-----------------------------------------------------
   % Heat flows to vapor and liquid phases in (ET)
   %----------------------------------------------------- 
   % heat flow to vapor phase in (ET)
   QdotV2 = QdotWV2 - QdotVS2 - pdV2 ...   %+ P.ratio_top_bottom*Jtr * X_v * h_X_v ... %(htr+0.5*vtr^2)... % term added for top fill
        + P.ratio_top_bottom * Jtr * (htr_L+0.5*vtr^2-qh2)... % term added for top fill
        - Jvvalve2*(hvalve2 + 0.5*vv2^2) ...
        - Jcd2*hcd2;
    
    % heat flow to liquid phase in (ET) % added by GP
    QdotL2 = QdotWL2 - QdotLS2 + pdV2 ...
           + (1-P.ratio_top_bottom)*Jtr*(htr_L+0.5*vtr^2) ...%+ P.ratio_top_bottom * Jtr* (1-X_v)  * h_X_L ...
           + Jcd2*hcd2;
  
    %-----------------------------------------------------
    % Variation of internal energies (boundary layers and bulk)
    %-----------------------------------------------------    
    % internal energies of liquid boundary layers in (ST)
    for i=1:P.nL1-1
        if i==1
            TL1im1 = Ts1;
        else
            TL1im1 = TL1(i-1);
        end
        rho_L1i = refpropm('D','T',TL1(i),'Q',0,'PARAHYD');
        duL1dt(i) = ((TL1(i+1)-TL1(i))/l12_L1(i+1)-(TL1(i)-TL1im1)/l12_L1(i))*P.kappa_L /(l_L1(i)*rho_L1i);
    end   
    % ENERGY BALANCE FOR THE LIQUID IN (ST)
    duL1dt(P.nL1) = (QdotL1 - JL1*(refpropm('U','T',TL1(P.nL1),'Q',0,'PARAHYD')))/mL1;
       
    %  internal energies of vapor boundary layers in  (ST)
    for i=1:P.nV1-1
        if i==1
            Tv1im1 = Ts1;
        else
            Tv1im1 = Tv1(i-1);
        end
        rhov1i= refpropm('D','P',pv1/1000,'U',uv1(i),'PARAHYD');
        duv1dt(i) = ((Tv1(i+1)-Tv1(i))/l12_V1(i+1)-(Tv1(i)-Tv1im1)/l12_V1(i))*P.kappa_v /(l_V1(i)*rhov1i);
    end
    % ENERGY BALANCE FOR THE VAPOR IN (ST)
    duv1dt(P.nV1) = (QdotV1 - Jv1*(refpropm('U','T',Tv1(P.nV1),'D',rhov1,'PARAHYD')))/mv1;    
    
    %  internal energies of liquid boundary layers in  (ET)
    for i=1:P.nL2-1
        if i==1
            TL2im1 = Ts2;
        else
            TL2im1 = TL2(i-1);
        end
        rho_L2i = refpropm('D','T',TL2(i),'Q',0,'PARAHYD');
        duL2dt(i) = ((TL2(i+1)-TL2(i))/l12_L2(i+1)-(TL2(i)-TL2im1)/l12_L2(i))*P.kappa_L /(l_L2(i)*rho_L2i);
    end
    % ENERGY BALANCE FOR THE BULK LIQUID IN (ET)
    duL2dt(P.nL2) = (QdotL2 - JL2*(refpropm('U','T',TL2(P.nL2),'Q',0,'PARAHYD')))/mL2;
    if mL2<=0.9 % Tuning for uL2 smoothness
        duL2dt(:) = 0;
    end
    
    
    %  internal energies of vapor boundary layers in  (ET)
    for i=1:P.nV2-1
        if i==1
            Tv2im1 = Ts2;
        else
            Tv2im1 = Tv2(i-1);
        end
        rhov2i= refpropm('D','P',pv2/1000,'U',uv2(i),'PARAHYD');
        duv2dt(i) = ((Tv2(i+1)-Tv2(i))/l12_V2(i+1)-(Tv2(i)-Tv2im1)/l12_V2(i))*P.kappa_v /(l_V2(i)*rhov2i);
    end
    % ENERGY BALANCE FOR THE BULK VAPOR IN (ET)
    duv2dt(P.nV2) = (QdotV2 - Jv2*(refpropm('U','T',Tv2(P.nV2),'D',rhov2,'PARAHYD')))/mv2;

    
    % ET wall temperature
    cw2 = 2.516173240451E-11*Tw2^6 - 2.695483209737E-08*Tw2^5 + 0.00001122596286143*Tw2^4 - 0.002261465800734*Tw2^3 + 0.214810433559*Tw2^2 - 5.41715155529*Tw2^1 + 51.75489930095; % temperature dependent specific heat capacity for stainless steel 304, in J/K
    dcw2dT=6*2.516173240451E-11*Tw2^5 -5* 2.695483209737E-08*Tw2^4 + 4*0.00001122596286143*Tw2^3 -3* 0.002261465800734*Tw2^2 + 2*0.214810433559*Tw2- 5.41715155529; % derivative of the specific heat as a function of temperature.
    P.QdotEW2 = -7.462776654302E-02*VL2^2 + 4.445867251697E+00*VL2 + 3.108170556297E+01;    % correlation for heat transfer profile of 3,300 gallon vertical Dewar at LLNL. VL2 in m^3.
    dTw2dt = (P.QdotEW2 - QdotWL2 - QdotWV2)/(P.mw2 * ( cw2 + Tw2 *dcw2dT));                % variation of wall temperature, including actual temperature dependent heat capacity
        
    % variables to be used for post-processing
    AAA = Jvvalve1;
    BBB = - QdotVS1;
    CCC = pdV1;
    DDD = Jvvalve1*(hvalve1+0.5*vv1^2);
    EEE =  - Jcd1*hcd1;
    FFF = Jboil*hboil;
    GGG = - QdotLS1;
    HHH = - Jtr*htr_L+0.5*vtr^2;
    III = + Jcd1*hcd1 ;
    JJJ = - Jvap*htr_L;
    KKK = QdotV1 ;
    LLL = QdotL1 ;
    
    MMM = QdotWV2;
    NNN = - QdotVS2;
    OOO = pdV2;   
    PPP =  P.ratio_top_bottom * Jtr * (htr_L+0.5*vtr^2-qh2);
    QQQ = - Jvvalve2*(hvalve2 + 0.5*vv2^2);
    RRR = - Jcd2*hcd2;
    SSS = QdotWL2;
    TTT = - QdotLS2;
    UUU = + (1-P.ratio_top_bottom)*Jtr*(htr_L+0.5*vtr^2);
    VVV =  + Jcd2*hcd2;
    WWW = QdotV2;
    XXX = QdotL2;
   
    % state derivatives
    dxdt(1) = JL1;
    dxdt(2:1+P.nL1) = duL1dt;
    dxdt(P.nL1+2) = Jv1;
    dxdt(P.nL1+3:P.nL1+2+P.nV1) =duv1dt;
    dxdt(P.nL1+P.nV1+3) = dTs1dt;
    dxdt(P.nL1+P.nV1+4) = dJtrdt;
    dxdt(P.nL1+P.nV1+5) = dmVapdt;
    dxdt(P.nL1+P.nV1+6) = dJboildt;
    dxdt(P.nL1+P.nV1+7) = JL2;
    dxdt(P.nL1+P.nV1+8:P.nL1+P.nV1+P.nL2+7) = duL2dt;
    dxdt(P.nL1+P.nV1+P.nL2+8) = Jv2;
    dxdt(P.nL1+P.nV1+P.nL2+9:P.nL1+P.nV1+P.nL2+P.nV2+8) = duv2dt;
    dxdt(P.nL1+P.nV1+P.nL2+P.nV2+9) = dTs2dt;
    dxdt(P.nL1+P.nV1+P.nL2+P.nV2+10) = dTw2dt;
    
    % the following is only so that some variables can be saved
    dxdt(P.nL1+P.nV1+P.nL2+P.nV2+11) = 1*(Jcd1-x(P.nL1+P.nV1+P.nL2+P.nV2+11));
    dxdt(P.nL1+P.nV1+P.nL2+P.nV2+12) = 1*(Jcd2-x(P.nL1+P.nV1+P.nL2+P.nV2+12));
    dxdt(P.nL1+P.nV1+P.nL2+P.nV2+13) = 1*(AAA-x(P.nL1+P.nV1+P.nL2+P.nV2+13));
    dxdt(P.nL1+P.nV1+P.nL2+P.nV2+14) = 1*(BBB-x(P.nL1+P.nV1+P.nL2+P.nV2+14));
    dxdt(P.nL1+P.nV1+P.nL2+P.nV2+15) = 1*(CCC-x(P.nL1+P.nV1+P.nL2+P.nV2+15));
    dxdt(P.nL1+P.nV1+P.nL2+P.nV2+16) = 1*(DDD-x(P.nL1+P.nV1+P.nL2+P.nV2+16));
    dxdt(P.nL1+P.nV1+P.nL2+P.nV2+17) = 1*(EEE-x(P.nL1+P.nV1+P.nL2+P.nV2+17));
    dxdt(P.nL1+P.nV1+P.nL2+P.nV2+18) = 1*(FFF-x(P.nL1+P.nV1+P.nL2+P.nV2+18));
    dxdt(P.nL1+P.nV1+P.nL2+P.nV2+19) = 1*(GGG-x(P.nL1+P.nV1+P.nL2+P.nV2+19));
    dxdt(P.nL1+P.nV1+P.nL2+P.nV2+20) = 1*(HHH-x(P.nL1+P.nV1+P.nL2+P.nV2+20));
    dxdt(P.nL1+P.nV1+P.nL2+P.nV2+21) = 1*(III-x(P.nL1+P.nV1+P.nL2+P.nV2+21));
    dxdt(P.nL1+P.nV1+P.nL2+P.nV2+22) = 1*(JJJ-x(P.nL1+P.nV1+P.nL2+P.nV2+22));
    dxdt(P.nL1+P.nV1+P.nL2+P.nV2+23) = 1*(KKK-x(P.nL1+P.nV1+P.nL2+P.nV2+23));
    dxdt(P.nL1+P.nV1+P.nL2+P.nV2+24) = 1*(LLL-x(P.nL1+P.nV1+P.nL2+P.nV2+24));
    dxdt(P.nL1+P.nV1+P.nL2+P.nV2+25) = 1*(MMM-x(P.nL1+P.nV1+P.nL2+P.nV2+25));
    dxdt(P.nL1+P.nV1+P.nL2+P.nV2+26) = 1*(NNN-x(P.nL1+P.nV1+P.nL2+P.nV2+26));
    dxdt(P.nL1+P.nV1+P.nL2+P.nV2+27) = 1*(OOO-x(P.nL1+P.nV1+P.nL2+P.nV2+27));
    dxdt(P.nL1+P.nV1+P.nL2+P.nV2+28) = 1*(PPP-x(P.nL1+P.nV1+P.nL2+P.nV2+28));
    dxdt(P.nL1+P.nV1+P.nL2+P.nV2+29) = 1*(QQQ-x(P.nL1+P.nV1+P.nL2+P.nV2+29));
    dxdt(P.nL1+P.nV1+P.nL2+P.nV2+30) = 1*(RRR-x(P.nL1+P.nV1+P.nL2+P.nV2+30));
    dxdt(P.nL1+P.nV1+P.nL2+P.nV2+31) = 1*(SSS-x(P.nL1+P.nV1+P.nL2+P.nV2+31));
    dxdt(P.nL1+P.nV1+P.nL2+P.nV2+32) = 1*(TTT-x(P.nL1+P.nV1+P.nL2+P.nV2+32));
    dxdt(P.nL1+P.nV1+P.nL2+P.nV2+33) = 1*(UUU-x(P.nL1+P.nV1+P.nL2+P.nV2+33));
    dxdt(P.nL1+P.nV1+P.nL2+P.nV2+34) = 1*(VVV-x(P.nL1+P.nV1+P.nL2+P.nV2+34));
    dxdt(P.nL1+P.nV1+P.nL2+P.nV2+35) = 1*(WWW-x(P.nL1+P.nV1+P.nL2+P.nV2+35));
    dxdt(P.nL1+P.nV1+P.nL2+P.nV2+36) = 1*(XXX-x(P.nL1+P.nV1+P.nL2+P.nV2+36));

    % must return a column vector
    dxdt = dxdt';

    % update model structure (to store vent valve state)
    assignin('base','LH2Model',P);
   
end

    while tout(end) < P.tFinal
    %Solve until the first terminal event
        refine = 4;
        VentEvent = @(t,x) VentEvents(x,P,ETTVentState);

        try
            nt = length(t);
            options = odeset('MaxStep',1,'RelTol',P.relTol,'Events',VentEvent,'OutputSel',1,'Refine',refine,'InitialStep',t(nt)-t(nt-refine),'MaxStep',t(nt)-t(1));
        catch
            options = odeset('MaxStep',1,'RelTol',P.relTol,'Events',VentEvent,'OutputSel',1,'Refine',refine);
        end
        
        rhs = @(t,x) LH2dxdt(P,t,x);
        [t,x,te,xe,ie]= ode15s(rhs,[tstart,tfinal],x0,options);
  
        % Accumulate output.  This could be passed out as output arguments.
        nt = length(t);
        tout = [tout; t(2:nt)];
        ventstate =  ETTVentState * ones(nt,1);
        x = horzcat(x, ventstate); % concatenate output results from the ODE to the vent-state on (ET), 0 or 1
        xout = [xout; x(2:nt,:)];  % concatenate output results with previous results on different time windows.
        teout = [teout; te];       % Events at tstart are never reported.
        xeout = [xeout; xe];
        ieout = [ieout; ie];

        x0=x(end,:)';
        x0=x0(1:end-1); % last column (ETTVentState) is removed
   
        tstart = t(nt);
        ETTVentState = abs(ETTVentState - 1);
    end
    
    % close waitbar
    close(P.waitbar);
    display('Done with ODE solver');
    
    % configure data struct
    data.name = name;
    data.t = tout;

    % extract state variables
    data.mL1 = xout(:,1);
    data.uL1 = xout(:,2:1+P.nL1);
    data.mv1 = xout(:,P.nL1+2);
    data.uv1 = xout(:,P.nL1+3:P.nL1+2+P.nV1);
    data.Ts1 = xout(:,P.nL1+P.nV1+3);
    data.Jtr = xout(:,P.nL1+P.nV1+4);
    data.mVap = xout(:,P.nL1+P.nV1+5);
    data.Jboil = xout(:,P.nL1+P.nV1+6);
    data.mL2 = xout(:,P.nL1+P.nV1+7);
    data.uL2 = xout(:,P.nL1+P.nV1+8:P.nL1+P.nV1+P.nL2+7);
    data.mv2 = xout(:,P.nL1+P.nV1+P.nL2+8);
    data.uv2 = xout(:,P.nL1+P.nV1+P.nL2+9:P.nL1+P.nV1+P.nL2+P.nV2+8);
    data.Ts2 = xout(:,P.nL1+P.nV1+P.nL2+P.nV2+9);
    data.Tw2 = xout(:,P.nL1+P.nV1+P.nL2+P.nV2+10);
    data.Jcd1 = xout(:,P.nL1+P.nV1+P.nL2+P.nV2+11);
    data.Jcd2 = xout(:,P.nL1+P.nV1+P.nL2+P.nV2+12);
    
    data.AAA = xout(:,P.nL1+P.nV1+P.nL2+P.nV2+13);
    data.BBB = xout(:,P.nL1+P.nV1+P.nL2+P.nV2+14);
    data.CCC = xout(:,P.nL1+P.nV1+P.nL2+P.nV2+15);
    data.DDD = xout(:,P.nL1+P.nV1+P.nL2+P.nV2+16);
    data.EEE = xout(:,P.nL1+P.nV1+P.nL2+P.nV2+17);
    data.FFF = xout(:,P.nL1+P.nV1+P.nL2+P.nV2+18);
    data.GGG = xout(:,P.nL1+P.nV1+P.nL2+P.nV2+19);
    data.HHH = xout(:,P.nL1+P.nV1+P.nL2+P.nV2+20);
    data.III = xout(:,P.nL1+P.nV1+P.nL2+P.nV2+21);
    data.JJJ = xout(:,P.nL1+P.nV1+P.nL2+P.nV2+22);
    data.KKK = xout(:,P.nL1+P.nV1+P.nL2+P.nV2+23);
    data.LLL = xout(:,P.nL1+P.nV1+P.nL2+P.nV2+24);
    data.MMM = xout(:,P.nL1+P.nV1+P.nL2+P.nV2+25);
    data.NNN = xout(:,P.nL1+P.nV1+P.nL2+P.nV2+26);
    data.OOO = xout(:,P.nL1+P.nV1+P.nL2+P.nV2+27);
    data.PPP = xout(:,P.nL1+P.nV1+P.nL2+P.nV2+28);
    data.QQQ = xout(:,P.nL1+P.nV1+P.nL2+P.nV2+29);
    data.RRR = xout(:,P.nL1+P.nV1+P.nL2+P.nV2+30);
    data.SSS = xout(:,P.nL1+P.nV1+P.nL2+P.nV2+31);
    data.TTT = xout(:,P.nL1+P.nV1+P.nL2+P.nV2+32);
    data.UUU = xout(:,P.nL1+P.nV1+P.nL2+P.nV2+33);
    data.VVV = xout(:,P.nL1+P.nV1+P.nL2+P.nV2+34);
    
    data.WWW = xout(:,P.nL1+P.nV1+P.nL2+P.nV2+35);
    data.XXX = xout(:,P.nL1+P.nV1+P.nL2+P.nV2+36);
    
    data.ETTTVenstate = xout(:,P.nL1+P.nV1+P.nL2+P.nV2+37); % ET Vent State
end

function y = dsqrt(x)
% directed square root
y = sqrt(abs(x)).*sign(x);
end

function mdot = gasFlow(CA,gamma,rho,P1,P2)
% choked/nonchoked flow. 
if P1<P2
	mdot = -gasFlow(CA,gamma,rho,P2,P1);
else
	%assumes P1 always >= P2
	threshold = ((gamma+1)/2)^(gamma/(gamma-1));
	if P1/P2 >= threshold
		% choked flow
		mdot = CA*sqrt(gamma*rho*P1*(2/(gamma+1))^((gamma+1)/(gamma-1)));
	else
		% nonchoked
		mdot = CA*sqrt(2*rho*P1*(gamma/(gamma-1))*((P2/P1)^(2/gamma)-(P2/P1)^((gamma+1)/gamma)));
	end
end
end


function H=cylVToH(V,R,L) 
%Cylinder height function
    A=pi*(R^2);s=V/L; x=0.01;error=1; b=1;
    if s>A/2
        sup=abs(s-A);
    else
        sup=s;
    end
    fun=@(x) sup-b*(R^2)*atan((((R^2)-(x^2))^(1/2))/x)+x*((R^2)-(x^2))^(1/2);
    
    while error>=1e-4;
        xold=x;  y=((R^2)-(x^2))^(1/2); f=fun(x);
        alpha=-(((R^2)-(x^2))^(-1/2))-(((R^2)-(x^2))^(1/2))/(x^2);
        supd=-b*alpha*(R^2)*(1/(((y/x)^2)+1))+(((R^2)-(x^2))^(1/2))-(x^2)*(((R^2)-(x^2))^(-1/2));
        x=x-f/supd; error=abs((xold-x)/xold);
    end
    if s>=A/2
        H=R+x;
    else
        H=R-x;
    end
end

function pv=vaporpressure(uv,rhov)
% new function for vapor pressure, based on internal energy and density
  try
       quality=refpropm('q','D',rhov,'U',uv,'PARAHYD');
  catch
      uv=fix(100*uv)/100;
      try
          quality=refpropm('q','D',rhov,'U',fix(100*uv)/100,'PARAHYD');
          display('Non-truncation of uv did not converge in "quality". Choosing truncated value instead');
      catch
          quality=refpropm('q','D',fix(100*rhov)/100,'U',fix(100*uv)/100,'PARAHYD');
          rhov=fix(100*rhov)/100;
        %uv=fix(100*uv)/100;
      end
  end
    if quality < 1 && quality >0
         % 2 phase
         temp=refpropm('T','D',rhov,'U',uv,'PARAHYD');
         pv = refpropm('P','T',temp,'Q',1,'PARAHYD')*1e3;% return value in Pa
    else
        %supercricical
        pv = refpropm('P','D',rhov,'U',uv,'PARAHYD')*1e3;% return value in Pa
    end
end


function [value,isterminal,direction] = VentEvents(x,P,ventstate) 
% Stops ODE solver every time the state of the vent valve in (ET) changes
rho_L22= -5.12074746E-07*(x(P.nL1+P.nV1+P.nL2+7)/1000)^3 - 1.56628367E-05*(x(P.nL1+P.nV1+P.nL2+7)/1000)^2 - 1.18436797E-01*(x(P.nL1+P.nV1+P.nL2+7)/1000) + 7.06218354E+01;
VL22 = x(P.nL1+P.nV1+7)/rho_L22;
Vullage22 = P.VTotal2-VL22;
rhov22 = x(P.nL1+P.nV1+P.nL2+8)/Vullage22;
p22=vaporpressure(x(P.nL1+P.nV1+P.nL2+P.nV2+8),rhov22);
value = [p22-P.p_ET_low; p22-P.p_ET_high]; % The value that we want to be zero

if ventstate > 0;          % this is to make sure the vent valve does not open when the pressure is dropping in ET
    isterminal = [1 ; 1];  % Halt integration if = 1. If venting is taking place, alway halt integration  
else
    isterminal = [0 ; 1];  % Halt integration if = 1. If no venting, then halt integration only if venting pressure is reached
end

direction =  [-1; +1];     % value of +1 locates only zeros where the event function is increasing, and -1 locates only zeros where the event function is decreasing.

end

