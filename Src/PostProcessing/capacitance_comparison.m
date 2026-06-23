% Codice per creare grafico con confronto capacita, con problemi di
% dimensioni diverse.

eps_0 = 8.85e-12;
eps_r = 1;
rapporto_LX_LY = [0.1,0.5,1,2,10,20]; %Lx/Ly
C_ideal = eps_0*eps_r*rapporto_LX_LY; % * Lx/Ly

C_charge = [6.610893e-12,1.179049e-11,1.711091e-11,2.736406e-11,1.043425e-10,1.977562e-10];
C_energy = [7.62572e-12,1.34677e-11,1.92484e-11,2.9187e-11,1.03379e-10,1.93202e-10];
figure;

%plot(rapporto_LX_LY, C_ideal, '-o', 'LineWidth', 2, 'MarkerSize', 7); hold on;
%plot(rapporto_LX_LY, C_charge, '-s', 'LineWidth', 2, 'MarkerSize', 7);
%plot(rapporto_LX_LY, C_energy, '-^', 'LineWidth', 2, 'MarkerSize', 7);

semilogy(rapporto_LX_LY, C_ideal, '-o', 'LineWidth', 2, 'MarkerSize', 7); 
hold on;
semilogy(rapporto_LX_LY, C_charge, '-s', 'LineWidth', 2, 'MarkerSize', 7);
semilogy(rapporto_LX_LY, C_energy, '-^', 'LineWidth', 2, 'MarkerSize', 7);

grid on;

xlabel('L_x/L_y');
ylabel('Capacitance (F/m)');
title('Capacitance vs L_x/L_y');

legend('C\_ideal','C\_charge','C\_energy','Location','best');

