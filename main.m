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
% Variables
Tin = 509; % not the official value, will be adjusted to achieve outlet of 533k 
Pin = 1000; %kPa
Tcin = 303;
D = .0245; %m
Lr = 15; % m
nflowinit = [.000009 .00041 .00005 0 0 .0000001 0 0]; % Same indices as pp
ntotthang = 0;
for i = 1:length(nflowinit)
    ntotthang = ntotthang + nflowinit(i);
end
molfracthang = nflowinit ./ntotthang;
disp(molfracthang)

% Set
Vr = D^2*pi/4*Lr; % m^3
numElements = 200;
MW = [0.02805, 0.03646, 0.01600, 0.1334, 0.04401, 0.0709, 0.09896, 0.01802]; %kg/mol
Htot = [136.77, 119.87, 1410.08, -89.92]; %kJ/mol
Cp = [0.0538, 0.0292, 0.0304, 0.1050, 0.0418, 0.0351, 0.0937, 0.0345]; % averaged over temp range, kJ/mol K
phi = .4;


%%%%%%%%%
% Logic %
%%%%%%%%%
disp(Vr)
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
mflowprod = ysoln(numElements,7)*MW(7);
disp("Kg flow product, 1 tube: " + num2str(mflowprod))
disp("N tubes: " + num2str(15.85/mflowprod));
plotdata(v, ysoln, conv);


