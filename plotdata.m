function null = plotdata(v, y, conv)
% Figure 1 -- flowrate vs reactor vol
figure(1)
f1 = y(:,1);

plot(v,y(:,1),v,y(:,2),v,y(:,3),v,y(:,4),v,y(:,5),'g-',v,y(:,6),v,y(:,7),v,y(:,8)) % might not work, check #color
grid
xlabel('Reactor Volume - L')
ylabel('Molar Flowrate - mol/hr')
title('Molar Flowrate vs. Reactor Volume')
legend('C_2H_4','HCl','O_2','C_2H_3Cl_3','CO_2','Cl_2','C_2H_4Cl_2','H2O','Location','northeastoutside')


% Figure 2 -- reactor T vs reactor vol
figure(2)
plot(v,y(:,9),'k-')
grid
xlabel('Reactor Volume - L')
ylabel('Reactor Temperature - K')
title('Reactor Temperature vs. Reactor Volume')

% Figure 3 -- coolant T vs reactor vol
figure(3)
plot(v,y(:,11),'k-')
grid
xlabel('Reactor Volume - L')
ylabel('Coolant Temperature - K')
title('Coolant Temperature vs. Reactor Volume')

% Figure 4 -- reactor P vs reactor vol
figure(4)
plot(v,y(:,10),'k-')
grid
xlabel('Reactor Volume - L')
ylabel('Reactor Pressure - kPa')
title('Reactor Pressure vs. Reactor Volume')

figure(5)
plot(v,conv,'k-')
grid
xlabel('Reactor Volume - L')
ylabel('Conversion (% of C_2H_4)')
title('Conversion profile')

end