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
nflowinit = [1,2,1,0,0,.125,0,0]; % Same indices as pp
phi = .4;
MW = [0.02805, 0.03646, 0.01600, 0.1334, 0.04401, 0.0709, 0.09896, 0.01802];%kg/mol
numElements = 200;

%%%%%%%%%
% Logic %
%%%%%%%%%
vspan = linspace(0, Vr, numElements);

%Loading Dependent variables
y0 = [nflowinit(1) nflowinit(2) nflowinit(3) nflowinit(4) nflowinit(5) nflowinit(6) nflowinit(7) nflowinit(8) Tin Pin];

handleranon = @(v,y) handler(v,y,phi,MW,Htot,Cp,(Vr/numElements));

%Mass Matrix (I have no idea if this is right, just assuming all diff eqs
%for now.)

M = [1     0     0     0     0     0     0     0     0     0;
     0     1     0     0     0     0     0     0     0     0;
     0     0     1     0     0     0     0     0     0     0;
     0     0     0     1     0     0     0     0     0     0;
     0     0     0     0     1     0     0     0     0     0;
     0     0     0     0     0     1     0     0     0     0;
     0     0     0     0     0     0     1     0     0     0;
     0     0     0     0     0     0     0     1     0     0;
     0     0     0     0     0     0     0     0     1     0;
     0     0     0     0     0     0     0     0     0     1;];
     

options = odeset('Mass',M);
%[ v, ysoln ] = ode45(handleranon,vspan,y0,options);
[ v, ysoln ] = ode45(handleranon,vspan,y0);

%plotdata();


