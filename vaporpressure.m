%Added by JMB
function pv=vaporpressure(Tv,rhov) % return pressure depending whether 2 phase or single phase
    if rhov<0;
        rhov=0.0001;
    end
    
    T_sat = -3.9389254667e-09*(rhov^6)+1.0053641879e-06*(rhov^5)...
        -1.0304184083e-04*(rhov^4)+5.3058942923e-03*(rhov^3)...
        -1.4792439609e-01*(rhov^2)+2.2234419496*rhov+1.7950995359e+01;
    % display(T_sat);
    % display(Tv);
     %display(rhov);
      %display('before if condition');
 if Tv>T_sat; % the fluid is in supercritical
     if Tv < 32.938 % this "if loop" is to avoid rounding error for refprop
         if Tv > 32.937
             Tv = 32.937;
      %       pause;
         end
     end 
     %display('Tv just before refprop function');
     %display(Tv);
     %display(rhov);
     %display('Refprop Calculation SC');
     %pv = (1.6133821043e-1*(Tv^3)-6.9432088540*Tv^2+1.1373052580e2*Tv-6.9558797798e2)*1e+3;
     pv = refpropm('P','T',Tv,'D',rhov,'PARAHYD')*1e3; % return value in Pa
    
 else % the fluid is in 2 phases
  %       display('Refprop Calculation 2 Phases');
     pv = (1.6133821043e-1*(Tv^3)-6.9432088540*Tv^2+1.1373052580e2*Tv-6.9558797798e2)*1e+3;
 end
 %display('end of if condition');
end