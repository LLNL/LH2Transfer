function pv=vaporpressure(uv,rhov)
% new function for vapor pressure, based on internal energy and density
  try
       quality=refpropm('q','D',rhov,'U',uv,'PARAHYD');
  catch
      uv=fix(100*uv)/100;
      try
          quality=refpropm('q','D',rhov,'U',fix(100*uv)/100,'PARAHYD');
          DISP('Non-truncation of uv did not converge in "quality". Choosing truncated value instead');
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
