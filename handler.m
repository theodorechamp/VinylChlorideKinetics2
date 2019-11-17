function outvar = handler(v, y, phi, MW, Htot, Cp, Lr, D)

T = y(9); %K
P_t = y(10); %kPa
Tc = y(11); %K

%Calculate total mole flow
nflowTot = 0; %mols
for i = 1:(length(y)-3)
   nflowTot = nflowTot + y(i);
end

%calculate mol fractions
molFracs = [0,0,0,0,0,0,0,0];
for i = 1:length(molFracs)
    molFracs(i) = y(i)/nflowTot;
end

%calculate partial pressures with Daltons law
pp = [0,0,0,0,0,0,0,0];
for i = 1:length(pp)
    pp(i) = molFracs(i) * P_t;
end

R = 8.3144621;

%Calculate rate constant for each reaction
a1 = 10^4.2; a2 = 10^13.23; a3 = 10^6.78;  % a's are the pre-exponential factors from the Lakshmanan paper.


% E's are acivation energies in units of J from the Lakshmanan paper.
E1 = -40100; E2 = -128080; E3 = -112000;

%Calculate rate equations.
k(1) = a1 * exp(E1/(R*T));
k(2) = a2 * exp(E2/(R*T));
k(3) = a3 * exp(E3/(R*T));
k(4) = (1000 * exp(17.13 - 13000/(1.987*T))) / exp(5.4+16000/(1.987*T));
%Rate constant units vary, see Lakshmanan paper for units

r1 = k(1) * pp(1) * pp(6)^0.5/3600;
r2 = k(2) * pp(7) * pp(6)^0.5/3600;
r3 = k(3) * pp(1) * pp(3) * pp(6)^0.5/3600;
r4 = k(4) * pp(3) / pp(6)/3600;
%mol/(L cat. * s)



outvar(1) = (-1*r1 - 1*r3) * (1-phi);
outvar(2) = (-2*r1 - 1*r2 + 4*r4) * (1-phi);
outvar(3) = (-0.5*r1 - 0.5*r2 - 3*r3 + r4) * (1-phi);
outvar(4) = r2 * (1-phi);
outvar(5) = (2*r3) * (1-phi);
outvar(6) = (2*r4) * (1-phi);
outvar(7) = (r1 - r2) * (1-phi);
outvar(8) = (r1 + r2 + 2*r3 + 2*r4) * (1-phi);

%summation of rate progression
term1 = r1*(1-phi)*Htot(1);
term1 = term1 + r2*(1-phi)*Htot(2);
term1 = term1 + r3*(1-phi)*Htot(3);
term1 = term1 + r4*(1-phi)*Htot(4);
term1 = term1 *1000; %kJ/(m^3*s) Conversion from L to m3

%Heat transfer parameters
U = .300; %kW/m2 k
term2 = U * 4 / D * (T - Tc); %kJ/(m^3*s)

%Heat transfer balance denominator
term3 = 0;
for i = 1:length(molFracs)
    term3 = term3 + y(i)*Cp(i); %Cp in kJ/mol K
end

outvar(9) = (term1 - term2) / term3; %K/m^3

nflow = [y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8)];
R = .0083144621; %m^3 kPa K-1 mol-1
rho = 0; % kg/m3 

for i = 1:length(nflow)
    rho = rho + pp(i)*MW(i)/(T*R); 
end

mflowTot = 0; % kg/s
for i = 1:length(MW)
    mflowTot = mflowTot + nflow(i)*MW(i);
end

Dp = D/8; %Using heuristic Weimer gave in class on 10/10/19 that particle diameter should be 1/8 of tube diameter
G0 = mflowTot/(pi*D^2/4); %kg/m^2 s
Ac = D^2 / 4 * pi; %m2
mu = 2.11*10^-5; %kg/(m*s) Values were averaged from Hysys simulation
outvar(10) = - G0/Dp*1/Ac*((1-phi)/phi^3)*((150*(1-phi)*mu)/Dp + 1.75*G0)*1/rho;
%outvar(10) = - 150/Ac*((1-phi)*mu/(Dp*G0)+(7/4))*((1-phi)/phi^3)*(G0^2/(rho*Dp))/1000; %kPa/m3 

Fc = 1; %kg/s
As = D * pi * Lr; %m2
Do = D + 2*.0036; %m
Cpc = ((Tc - 273) * .0029 + 1.5041 + 273)/1000; %kJ/kg K

outvar(11) = 4 * U * (T - Tc) / (Fc * Cpc * Do); %Our old method
%outvar(11) = pi()*Dp*mu*(Tc - T) / (Cpc * 

outvar = outvar';

end
