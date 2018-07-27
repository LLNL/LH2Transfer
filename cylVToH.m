    %Cylinder height function added by JMB
function H=cylVToH(V,R,L)
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