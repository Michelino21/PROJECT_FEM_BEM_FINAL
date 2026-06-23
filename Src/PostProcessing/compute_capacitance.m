function [C_energy_all,C_energy_inside,C_energy_clean] = compute_capacitance(V1,V2,e0,er,d_elm,conn,x,y,u,E_threshold,area,C_theoretical)
    
    %epsilon
    M = size(conn,1);
    eps_e = ones(M,1)*e0; eps_e(d_elm) = er*e0;

    %Electric Field
    Ex = zeros(M,1); Ey = zeros(M,1);
    for e = 1:M
        [Ex(e), Ey(e)] = Efield(e, x, y, conn, u);
    end
    Emag = sqrt(Ex.^2 + Ey.^2);

    %Good Element
    good_elm = d_elm(Emag(d_elm) <= E_threshold);

    C_energy_all = (1/(V1-V2)^2) * sum(eps_e .* (Ex.^2 + Ey.^2) .* area);

    C_energy_inside = (1/(V1-V2)^2) * sum(eps_e(d_elm) .* (Ex(d_elm).^2 + Ey(d_elm).^2) .* area(d_elm));

    C_energy_clean = (1/(V1-V2)^2) * sum(eps_e(good_elm) .* ...
               (Ex(good_elm).^2 + Ey(good_elm).^2) .* area(good_elm)); %solo elementi dentro e puliti.

    % Errore percentuale
    error_percent_energy_all = abs(C_energy_all - C_theoretical) / C_theoretical * 100;
    error_percent_energy_inside = abs(C_energy_inside - C_theoretical) / C_theoretical * 100;
    error_percent_energy_clean = abs(C_energy_clean - C_theoretical) / C_theoretical * 100;
    
    
    fprintf('\n');
    fprintf('------- CAPACITANCE ENERGY METHOD -------\n');
    fprintf('Capacitance in all the BEM domain       : %.6g F/m\n', C_energy_all);
    %fprintf('Capacitance inside the capacitor        : %.6g F/m\n', C_energy_inside);
    %fprintf('Capacitance inside the capacitor clean  : %.6g F/m\n', C_energy_clean);
    fprintf('Theoretical capacitance (infinite plate): %.6g F/m\n', C_theoretical);
    fprintf('Relative error all BEM                  : %.4g%%\n', error_percent_energy_all);
    %fprintf('Relative error inside capacitor         : %.4g%%\n', error_percent_energy_inside);
    %fprintf('Relative error inside capacitor         : %.4g%%\n', error_percent_energy_clean);
  
end