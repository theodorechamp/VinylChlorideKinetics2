% This file runs the script, ensuring the correct sequence of calling
% functions

%Species indices key:
    % 1 = c2h4
    % 2 = hcl
    % 3 = o2
    % 4 = 1,1,2-trichloroethane
    % 5 = co2
    % 6 = cl2
    % 7 = 1,2-dichloroethane
    % 8 = h2o

%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%
Tin = 533.579; % not the official value, will be adjusted to achieve outlet of 533 
Pin = 1000; %kPa
Tcin = 303;
Htot = [136.77, 119.87, 1410.08, -89.92]; %kJ/mol
Cp = [23.2, 15.8, 39.4, 15.5, 17.4, 19.7, 42.8, 18.0]; % averaged over temp range, kJ/mol K
D = .08; %m
Lr = 8; % m
Vr = D^2*pi/4*Lr; % m^3
nflowinit = [.00065,.0019,.000665,0,0,0.000001,0,0]; % Same indices as pp
phi = .4;
MW = [0.02805, 0.03646, 0.01600, 0.1334, 0.04401, 0.0709, 0.09896, 0.01802]; %kg/mol
numElements = 200;

%%%%%%%%%
% Logic %
%%%%%%%%%
vspan = linspace(0, Vr, numElements);

%Loading Dependent variables
y0 = [nflowinit(1) nflowinit(2) nflowinit(3) nflowinit(4) nflowinit(5) nflowinit(6) nflowinit(7) nflowinit(8) Tin Pin Tcin];

handleranon = @(v,y) handler(v,y,phi,MW,Htot,Cp,(Vr/numElements),Lr,D);

%for now.)     %Mass Matrix (I have no idea if this is right, just assuming all diff eqs


%options = odeset('Mass',M);
%[ v, ysoln ] = ode45(handleranon,vspan,y0,options);
[ v, ysoln ] = ode15s(handleranon,vspan,y0);
conv = zeros(numElements);
for i = 1:numElements
    conv(i) = (1-ysoln(i,1)/ysoln(1,1));
end
disp("Final Conversion: "+ num2str(conv(numElements)))
mflowprod = ysoln(200,7)*MW(7);
disp("Kg flow product, 1 tube: " + num2str(mflowprod))
disp("N tubes: " + num2str(15.85/mflowprod));
plotdata(v, ysoln, conv);


