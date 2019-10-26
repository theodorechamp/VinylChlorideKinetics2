% This file runs the script, ensuring the correct sequence of calling
% functions

%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%
Tin = 523; % not the official value, will be adjusted to achieve outlet of 533 
Pin = 1000; %kPa
Htot = [-136.77, -119.87, -1410.08, 89.92]; %kJ/mol
Cp = [2.32e4, 1.58e4, 3.94e4, 1.55e4, 1.74e4, 1.97e4, 4.28e4, 1.80e4]; % averaged over temp range
Vr = 1.414; % m^3
Lr = 3; % m
nflowinit = [10,20,10,0,0,1.25,0,0]; % Same indices as pp
phi = .4;
MW = [0.02805, 0.03646, 0.01600, 0.1334, 0.04401, 0.0709, 0.09896, 0.01802];%kg/mol
numElements = 200;

%%%%%%%%%
% Logic %
%%%%%%%%%
vspan = linspace(0, Vr, numElements);

y = zeros(numElements,10);

%Loading Dependent variables
y(1,1:8) = nflowinit;
y(1,9:10) = [Tin,Pin];


ysoln = euler1(y, phi, MW, Htot, Cp, (Vr/numElements), Vr);
%plotdata();


