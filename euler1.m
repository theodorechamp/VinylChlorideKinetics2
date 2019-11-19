function outvar = euler1(y,phi, MW, Htot, Cp, dV, Vr, D, Lr)

count = 1;

v = 0;
while v < Vr
    nflow = [y(count,1),y(count,2),y(count,3),y(count,4),y(count,5),y(count,6),y(count,7),y(count,8)];
    T = y(count,9);
    Pt = y(count,10);
    Tc = y(count,11);
    
    nflowTot = 0;
    for i = 1:8
       nflowTot = nflowTot + y(count,i);
    end
    
    molFracs = [0,0,0,0,0,0,0,0];
    for i = 1:length(molFracs)
        molFracs(i) = y(count,i)/nflowTot;
    end
    
    pp = [0,0,0,0,0,0,0,0];
    for i = 1:length(pp)
        pp(i) = molFracs(i) * Pt;
    end

    %init dydv with 0
    dydv = [0,0,0,0,0,0,0,0,0,0,0];

    RjmolK = 8.3144621; %J/mol K

    %Calculate rate constant for each reaction
    a1 = 10^4.2; a2 = 10^13.23; a3 = 10^6.78;  % a's are the pre-exponential factors from the Lakshmanan paper.


    % E's are acivation energies in units of J from the Lakshmanan paper.
    E1 = -40100; E2 = -128080; E3 = -112000;

    %Calculate rate equations.
    k(1) = a1 * exp(E1/(RjmolK*T));
    k(2) = a2 * exp(E2/(RjmolK*T));
    k(3) = a3 * exp(E3/(RjmolK*T));
    k(4) = (1000 * exp(17.13 - 13000/(1.987*T))) / exp(5.4+16000/(1.987*T));
    %Rate constant units vary, see Lakshmanan paper for units

    r1 = k(1) * pp(1) * pp(6)^0.5/3600;
    r2 = k(2) * pp(7) * pp(6)^0.5/3600;
    r3 = k(3) * pp(1) * pp(3) * pp(6)^0.5/3600;
    r4 = k(4) * pp(3) / pp(6)/3600;
    %mol/(L cat. * s)

    dydv(1) = (-1*r1 - 1*r3) * (1-phi);
    dydv(2) = (-2*r1 - 1*r2 + 4*r4) * (1-phi);
    dydv(3) = (-0.5*r1 - 0.5*r2 - 3*r3 + r4) * (1-phi);
    dydv(4) = r2 * (1-phi);
    dydv(5) = (2*r3) * (1-phi);
    dydv(6) = (2*r4) * (1-phi);
    dydv(7) = (r1 - r2) * (1-phi);
    dydv(8) = (r1 + r2 + 2*r3 + 2*r4) * (1-phi);

    
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

    dydv(9) = (term1 - term2) / term3; %K/m^3
    
    Rm3kpaKmol = .0083144621; %m^3 kPa K-1 mol-1
    
    rho = 0; % kg/m3 
    for i = 1:length(nflow)
        rho = rho + pp(i)*MW(i)/(T*Rm3kpaKmol); 
    end
    
    mflowTot = 0; % kg/s
    for i = 1:length(MW)
        mflowTot = mflowTot + nflow(i)*MW(i);
    end

    Dp = D/8; %Using heuristic Weimer gave in class on 10/10/19 that particle diameter should be 1/8 of tube diameter
    G0 = mflowTot/(pi*D^2/4); %kg/m^2 s
    Ac = D^2 / 4 * pi; %m2
    mu = 2.11*10^-5; %kg/(m*s) Values were averaged from Hysys simulation
    beta0 = (G0/rho/Dp) * ((1-phi)/phi^3) * (150*(1-phi)*mu/Dp+1.75*G0/1000); %kPa/m
    dydv(10) = - G0/Dp*1/Ac*((1-phi)/phi^3)*((150*(1-phi)*mu)/Dp + 1.75*G0)*1/rho;
    
    Fc = 1; %kg/s
    As = D * pi * Lr; %m2
    Do = D + 2*.0036; %m
    Cpc = ((Tc - 273) * .0029 + 1.5041 + 273)/1000; %kJ/kg K

    dydv(11) = 4 * U * (T - Tc) / (Fc * Cpc * Do);
    
    for i = 1:length(y(count,:))
       y(count+1,i) = y(count,i) + dydv(i)*dV;
    end
    
    v = v + dV;
    count = count + 1;
end

outvar = y;
end
