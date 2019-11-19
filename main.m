% This file runs the script, and allows values to be inputted
close all
clear
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
Tin = 471; % K 
Pin = 1000; %kPa
Tcin = 303; %K
D = .05; % m Using value from working solution
Lr = 12; % m
nflowinit = [.147 .432 .154 0 0 1.54e-8 0 0];
%mols/sec

% Set
Vr = D^2*pi/4*Lr; % m^3
numElements = 20; % number of solver iterations
dV = Vr/numElements;
MW = [0.02805, 0.03646, 0.01600, 0.1334, 0.04401, 0.0709, 0.09896, 0.01802]; %kg/mol
Htot = [-240.14, -149.31, -1321.08, -116.82]; %kJ/mol
Cp = [0.0538, 0.0292, 0.0304, 0.1050, 0.0418, 0.0351, 0.0937, 0.0345]; % averaged over temp range, kJ/mol K
phi = .4;


%%%%%%%%%
% Logic %
%%%%%%%%% 
y = zeros(numElements, 11);
vspan = linspace(0,Vr,numElements);

%Loading Dependent variables
y(1,:) = [nflowinit(1) nflowinit(2) nflowinit(3) nflowinit(4) nflowinit(5) nflowinit(6) nflowinit(7) nflowinit(8) Tin Pin Tcin];
ysoln = euler1(y,phi, MW, Htot, Cp, dV, Vr, D, Lr);


conv = zeros(numElements);
for i = 1:numElements
    conv(i) = (1-ysoln(i,1)/ysoln(1,1));
end

disp("Final Conversion: "+ num2str(conv(numElements)))
plotdata(vspan, ysoln, conv);


