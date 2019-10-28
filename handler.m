function outvar = handler(v, y, phi, MW, Htot, Cp, dV, Lr, D)

%unpacking y
%disp y(1)
%{
nflowc2h4 = y(1);
nflowhcl = y(2);
nflowo2 = y(3);
nflowc2h3cl3 = y(4);
nflowco2 = y(5);
nflowcl2 = y(6);
nflowc2h4cl2 = y(7);
nflowh2o = y(8);
%}

T = y(9); %K
P_t = y(10); %kPA
Tc = y(11); %K


%Check with Matt if mol flowrates or just mol comps should be used.
%calculate total mole flow
nflowTot = 0; %mols
for i = 1:(length(y)-3)
   nflowTot = nflowTot + y(i);
end
%disp(int2str(nflowTot) + "  = 300");

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

%a2 is very large, causes a very stiff system (very fast rate expression
% E's are acivation energies in units of J from the Lakshmanan paper.
E1 = -40100; E2 = -128080; E3 = -112000;

%Calculate rate equations.
k(1) = a1 * exp(E1/(R*T));
k(2) = a2 * exp(E2/(R*T));
k(3) = a3 * exp(E3/(R*T));
k(4) = (1000 * exp(17.13 - 13000/(1.987*T))) / exp(5.4+16000/(1.987*T));

r1 = k(1) * pp(1) * pp(6)^0.5/3600;
r2 = k(2) * pp(7) * pp(6)^0.5/3600;
r3 = k(3) * pp(1) * pp(3) * pp(6)^0.5/3600;
r4 = k(4) * pp(3) / pp(6)/3600;
%Rates units vary, now s^-1, but good because k and pp are good


outvar(1) = (-1*r1 - 1*r3) * (1-phi);
outvar(2) = (-2*r1 - 1*r2 + 4*r4) * (1-phi);
outvar(3) = (-0.5*r1 - 0.5*r2 - 3*r3 + r4) * (1-phi);
outvar(4) = r2 * (1-phi);
outvar(5) = (2*r3) * (1-phi);
outvar(6) = (2*r4) * (1-phi);
outvar(7) = (r1 - r2) * (1-phi);
outvar(8) = (r1 + r2 + 2*r3 + 2*r4) * (1-phi);
%check to make sure correct

%summation of rate progression
term1 = r1*(1-phi)*Htot(1);
term1 = term1 + r2*(1-phi)*Htot(2);
term1 = term1 + r3*(1-phi)*Htot(3);
term1 = term1 + r4*(1-phi)*Htot(4);


%Params
U = .300; %kW/m2 k

term2 = U * dV * 4 / D * (T - Tc);

term3 = 0;
for i = 1:length(molFracs)
    term3 = term3 + y(i)*Cp(i); %Cp in kJ/mol K
end

outvar(9) = (term1 - term2) / term3;

nflow = [y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8)];
R = .0083144621; %m3 kPa K-1 mol-1
rho = 0; % kg/m3

for i = 1:length(nflow)
    rho = rho + pp(i)*MW(i)/(T*R); 
end
%disp(rho)
mflowTot = 0; % kg/hr
for i = 1:length(MW)
    mflowTot = mflowTot + nflow(i)*MW(i);
end

%m
Dp = D/8; %Using heuristic Weimer gave in class on 10/10/19 that particle diameter should be 1/8 of tube diameter
G0 = mflowTot/(pi*D^2/4); %kg/m^2 s
Ac = D^2 / 4 * pi; %m2
mu = 2.11*10^-5; %From Hysys
outvar(10) = - 150/Ac*((1-phi)*mu/(Dp*G0)+(7/4))*((1-phi)/phi^3)*(G0^2/(rho*Dp))/1000; %kPa/m3 
%disp(G0*Ac/rho)
Fc = 3.15; %kg/s
As = D * pi * Lr;
Do = D + 2*.0036;
Cpc = ((Tc - 273) * .0029 + 1.5041 + 273)/1000; %kJ/kg K

%disp(Cpc)
outvar(11) = 4 * U * (T - Tc) / (Fc * Cpc * Do);

outvar = outvar';

end
