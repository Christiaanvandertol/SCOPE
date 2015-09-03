function fVcmo = calcfVcmo(P,x)
fVcmo = P(1)+(P(2)*x.^P(6)).*(1./(P(4)+P(3)*x.^P(5)));
