function outvar = handler(v, y, phi, MW, Htot, Cp, dV)
%testing branch
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

T = y(9);
P_t = y(10);


%Check with Matt if mol flowrates or just mol comps should be used.
%calculate total mole flow
nflowTot = 0;
for i = 1:(length(y)-2)
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
% E's are acivation energies in units of kJ from the Lakshmanan paper.
E1 = -40.1; E2 = -128.08; E3 = -112;

%Calculate rate equations.
k(1) = a1 * exp(E1/(R*T));
k(2) = a2 * exp(E2/(R*T));
k(3) = a3 * exp(E3/(R*T));
k(4) = (1000 * exp(17.13 - 13000/(1.987*T))) / exp(5.4+16000/(1.987*T));

r1 = k(1) * pp(1) * pp(6)^0.5;
r2 = k(2) * pp(7) * pp(6)^0.5;
r3 = k(3) * pp(1) * pp(3) * pp(6)^0.5;
r4 = k(4) * pp(3) / pp(6);

%summation of rate progression
term1 = r1*Htot(1);
term1 = term1 + r2*Htot(2);
term1 = term1 + r3*Htot(3);
term1 = term1 + r4*Htot(4);


%Params
U = 300; %Given in problem statement.
D = .0245;
Tinf = 298; %Dummy Value, 25C

term2 = U * dV * 4 / D * (y(9) - Tinf);

term3 = 0;
for i = 1:length(molFracs)
    term3 = term3 + molFracs(i)*Cp(i);
end

outvar(1) = (-1*r1 - 1*r3) * (1-phi);
outvar(2) = (-2*r1 - 1*r2 + 4*r4) * (1-phi);
outvar(3) = (-0.5*r2 - 3*r3 + r4) * (1-phi);
outvar(4) = r2 * (1-phi);
outvar(5) = (2*r3) * (1-phi);
outvar(6) = (2*r4) * (1-phi);
outvar(7) = (r1 - 1*r2) * (1-phi);
outvar(8) = (r1 + r2 + 2*r3 + 2*r4) * (1-phi);

outvar(9) = (term1 - term2) / term3;

nflow = [y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8)];
R = 8.3144621; %L kPa K-1 mol-1
rho = 0;

for i = 1:length(nflow)
    rho = rho + pp(i)*MW(i)/(T*R);
end

mflowTot = 0; % mol/hr
for i = 1:length(MW)
    mflowTot = mflowTot+nflow(i)*MW(i);
end

D = 0.0254;
Dp = D/8; %Using heuristic Weimer gave in class on 10/10/19 that particle diameter should be 1/8 of tube diameter
G0 = mflowTot/(pi*D^2/4);
Ac = D * pi;
mu = 2.11*10^-5; %From Hysys
outvar(10) = - 150/Ac*((1-phi)/(Dp*G0/(mu))+(7/4))*((1-phi)/phi^3)*(G0^2/(rho*Dp));

outvar = outvar';

end
