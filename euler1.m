function outvar = euler1(y,phi, MW, Htot, Cp, dV, vr)

%y(*,1) refers to y0
count = 1;

v = 0;
while v < vr
    T = y(count,9);
    Pt = y(count,10);
    
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
    dydv = [0,0,0,0,0,0,0,0,0,0];
    
    rates = rateODE(pp, T);
    dydv(1:8) = speciesBal(rates, phi);
    
    dydv(9) = energyBal(Htot, rates, T, molFracs, Cp, dV);
    dydv(10) = ergun(T, pp, MW, y(count,:), phi);
    
    
    for i = 1:length(y(count))
       y(count+1,i) = y(count,i) + dydv(i)*dV;
    end
    
    v = v + dV;
    count = count + 1;
end

outvar = y;
end
