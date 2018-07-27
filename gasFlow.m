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