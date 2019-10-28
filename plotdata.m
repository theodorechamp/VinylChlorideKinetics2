function null = plotdata(v, y, conv)
% Figure 1 -- flowrate vs reactor vol
figure(1)
%yyaxis left
plot(v.*1000,y(:,1),'b-',v.*1000,y(:,8),'y-',v.*1000,y(:,3),'m-',v.*1000,y(:,4),'c-',v.*1000,y(:,5),'r-',v.*1000,y(:,6),'g-',v.*1000,y(:,7),'k-')
ylabel('Molar Flowrate - mol/hr')
%yyaxis right
grid
xlabel('Reactor Volume - L')
ylabel('HCl Molar Flowrate - mol/hr')
title('Flowrate vs. Reactor Volume')
legend('C_2H_4','H_2O','O_2','C_2H_3Cl_3','CO_2','Cl_2','C_2H_4Cl_2','Location','northeastoutside')

figure(2)
plot(v.*1000,y(:,2),'k-')
grid
xlabel('Reactor Volume - L')
ylabel('HCl Molar Flowrate - mol/hr')
title('HCl Flowrate vs. Reactor Volume')


% Figure 3 -- reactor T vs reactor vol
figure(3)
plot(v.*1000,y(:,9),'k-')
grid
xlabel('Reactor Volume - L')
ylabel('Reactor Temperature - K')
title('Reactor Temperature vs. Reactor Volume')

% Figure 4 -- coolant T vs reactor vol
figure(4)
plot(v.*1000,y(:,11),'k-')
grid
xlabel('Reactor Volume - L')
ylabel('Coolant Temperature - K')
title('Coolant Temperature vs. Reactor Volume')

% Figure 5 -- reactor P vs reactor vol
figure(5)
plot(v.*1000,y(:,10),'k-')
grid
xlabel('Reactor Volume - L')
ylabel('Reactor Pressure - kPa')
title('Reactor Pressure vs. Reactor Volume')

figure(6)
plot(v.*1000,conv,'k-')
grid
xlabel('Reactor Volume - L')
ylabel('Conversion (% of C_2H_4)')
title('Conversion profile')

end