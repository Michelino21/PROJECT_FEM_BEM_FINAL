function [] = error_FEM(Lx,Ly,delx,dely,e0,er,V1,V2,u_FEM_BEM,t_FEM_BEM,folder_name,SALVATAGGIO,margin_BEM_coeff,margini_list,x_r, y_r,outb_r, ptop_r, pbottom_r,TIME_mesh_bem_generation,E_threshold,C_theoretical,C_energy_all,C_energy_inside,C_energy_clean)
    
    fprintf('\n==========================================\n');
    fprintf('               ERROR ANALISYS               \n');
    fprintf('==========================================\n');
   
    %% === STEP 1: SETUP
    n_test  = length(margini_list);
    
    err_V  = zeros(n_test, 1);
    err_C  = zeros(n_test, 1);
    err_C_energy_all   = zeros(n_test, 1);
    err_C_energy_inside   = zeros(n_test, 1);
    err_C_energy_clean   = zeros(n_test, 1);

    t_fem_mesh = zeros(n_test, 1);
    t_fem_solve = zeros(n_test, 1);
    
    N_nodes = zeros(n_test, 1);

    % creo funzione interpolante BEM, la interpolo su dei punti di test
    x_test = x_r;
    y_test = y_r;

    F_interp_bem = scatteredInterpolant(x_r, y_r, u_FEM_BEM, 'linear', 'nearest');
    u_bem_on_ref = F_interp_bem(x_test, y_test);
    
    %% === STEP 2: FEM PURO CON MARGINI CRESCENTI ===
    
    for k = 1:n_test
        fprintf('\n--- iteration: %d , margin: %d * Ly---\n',k,margini_list(k));
        
        % Dati FEM

        margine =  margini_list(k) * Ly;
        X_mesh_k = Lx + 2*margine;
        Y_mesh_k = Ly + 2*margine;

        fprintf('Real Margin: %.4f m\n',margine);
        fprintf('X_mesh           : %.4f m\n', X_mesh_k);
        fprintf('Y_mesh           : %.4f m\n', Y_mesh_k);
        
        % Genera mesh e risolvi FEM puro
        
        tic; 
        
        [conn_k, x_k, y_k, ~, ~, outb_k, ptop_k, pbottom_k, d_elm_k] = ...
            meshgen_cap(X_mesh_k, Y_mesh_k, Lx, Ly, delx, dely);
        
        t_fem_mesh(k) = toc;

        tic;
        
        [u_fem_k,area_k] = solve_FEM(conn_k, x_k, y_k, outb_k, ptop_k, pbottom_k, ...
                                  d_elm_k, e0, er, V1, V2);
        
        t_fem_solve(k) = toc;
        
        N_nodes(k) = length(x_k);

        % --- Plot potenziale FEM (Dominio Intero) (mappa colori) ---
        xu_fem = unique(x_k);
        yu_fem = unique(y_k);
        U_fem_res  = reshape(u_fem_k, length(xu_fem), length(yu_fem))';

        figure('Color', 'w');  % crea la figura con sfondo bianco
        clf;                   % pulisce la figura
        imagesc(xu_fem, yu_fem, U_fem_res);
        colormap(jet); axis equal tight;
        set(gca,'Ydir','normal'); colorbar;
        xlabel('x (m)'); ylabel('y (m)');
        title('Potential (V) - Pure FEM');
        subtitle(sprintf('FEM margin: %d * Ly', margini_list(k)));
        set(gcf,'Color',[1 1 1]); hold on
        plot(x_k(ptop_k),    y_k(ptop_k),    'k', 'linewidth', 5)
        plot(x_k(pbottom_k), y_k(pbottom_k), 'k', 'linewidth', 5)

        fig_name = sprintf('FEM_%d', margini_list(k));
        salva_fig(fig_name, folder_name, SALVATAGGIO);


        % Differenza FEM BEM
        
        % Interpola la soluzione FEM sul riferimento
        F_interp = scatteredInterpolant(x_k, y_k, u_fem_k, 'linear', 'nearest');
        u_fem_on_ref = F_interp(x_test, y_test);


        % --- Plot potenziale FEM on BEM domain ---
        xu = unique(x_test);
        yu = unique(y_test);
        U_fem_res_bem  = reshape(u_fem_on_ref, length(xu), length(yu))';

        figure('Color', 'w');  % crea la figura con sfondo bianco
        clf;                   % pulisce la figura
        imagesc(xu, yu, U_fem_res_bem);
        colormap(jet); axis equal tight;
        set(gca,'Ydir','normal'); colorbar;
        xlabel('x (m)'); ylabel('y (m)');
        title('Potential (V) - Pure FEM on BEM domain');
        subtitle(sprintf('FEM margin: %d * Ly', margini_list(k)));
        set(gcf,'Color',[1 1 1]); hold on
        plot(x_k(ptop_k),    y_k(ptop_k),    'k', 'linewidth', 5)
        plot(x_k(pbottom_k), y_k(pbottom_k), 'k', 'linewidth', 5)

        fig_name = sprintf('FEM_%d_BEM_%d', margini_list(k), margin_BEM_coeff);
        salva_fig(fig_name, folder_name, SALVATAGGIO);
        
        
        % Nodi interni (escludi bordo e piastre)
        inner = setdiff(1:length(x_r), [outb_r; ptop_r; pbottom_r]);

        % errore puntuale tra FEM e BEM (potenziale)
        err_u = u_fem_on_ref(inner) - u_bem_on_ref(inner);

        % Errore L2
        err_V(k) = sqrt(sum(err_u.^2) / sum(u_fem_on_ref(inner).^2));


        % Grafico errore puntuale
        err_u_on_ref = u_fem_on_ref - u_bem_on_ref;
        err_u_plot  = reshape(err_u_on_ref, length(xu), length(yu))';

        figure('Color', 'w');  % crea la figura con sfondo bianco
        clf;                    % pulisce la figura
        imagesc(xu, yu, abs(err_u_plot));
        colormap(jet); axis equal tight;
        set(gca,'Ydir','normal'); colorbar;
        xlabel('x (m)'); ylabel('y (m)');
        title(sprintf('Difference between u_{bem} and u_{fem} on BEM domain'));
        subtitle(sprintf('FEM margin: %d * Ly  |  BEM margin: %d * Ly', margini_list(k), margin_BEM_coeff));
        set(gcf,'Color',[1 1 1]); hold on
        plot(x_r(ptop_r),    y_r(ptop_r),    'k', 'linewidth', 5)
        plot(x_r(pbottom_r), y_r(pbottom_r), 'k', 'linewidth', 5)

        fig_name = sprintf('u_err_FEM_%d_BEM_%d', margini_list(k), margin_BEM_coeff);
        salva_fig(fig_name, folder_name, SALVATAGGIO);

        if(SALVATAGGIO == 1)
            %Salva le variabili ogni iterazione, per garantire salvataggio
            %parziale
            save(fullfile(folder_name, 'workspace_error.mat'));
        end
      
        % Errore sulla capacitanza

        [C_energy_all_k,C_energy_inside_k,C_energy_clean_k] = compute_capacitance(V1,V2,e0,er,d_elm_k,conn_k,x_k,y_k,u_fem_k,E_threshold,area_k,C_theoretical);

        err_C_energy_all(k) = C_energy_all_k;
        err_C_energy_inside(k) = C_energy_inside_k;
        err_C_energy_clean(k) = C_energy_clean_k;

        err_C(k) = abs(C_energy_all-C_energy_all_k)/C_energy_all_k*100;


        fprintf('Relative error all FEM-BEM vs FEM       : %.4g%%\n', err_C(k));

        if(SALVATAGGIO == 1)
            %Salva le variabili ogni iterazione, per garantire salvataggio
            %parziale
            save(fullfile(folder_name, 'workspace_error.mat'));
        end
        
        
    end
    
    %% === STEP 3: PLOT RISULTATI ===
    
    margini_label = margini_list;
    
    % --- Errore L2 vs margine ---
    figure('Color','w');
    semilogy(margini_label, err_V, 'bd-', 'LineWidth', 2, 'MarkerSize', 8)
    xlabel('margin / Ly');
    ylabel('Relative L^2 error on u');
    title('Pure FEM accuracy vs domain size (reference: FEM-BEM)');
    grid on
    
    salva_fig('L2_norm', folder_name, SALVATAGGIO);
    
    % --- Errore capacitanza metodo energia ---
    figure('Color','w');
    semilogy(margini_label, err_C, 'rs-', 'LineWidth', 2, 'MarkerSize', 8)
    xlabel('Domain margin / Ly');
    ylabel('Relative error on C (%), C_fem_bem/C_fem');
    title('Capacitance error energy method vs domain size');
    grid on
    
    salva_fig('C_energy', folder_name, SALVATAGGIO);

    % --- Tempo vs N nodi ---
    figure('Color','w');
    hold on
    %plot(margini_list, t_fem_mesh, 'gd-', 'LineWidth', 2, 'MarkerSize', 8, ...
    %     'DisplayName', 'Pure FEM mesh')
    %plot(margini_list, t_fem_solve, 'bd-', 'LineWidth', 2, 'MarkerSize', 8, ...
    %     'DisplayName', 'Pure FEM solve')
    plot(margini_list, t_fem_solve+t_fem_mesh, 'kd-', 'LineWidth', 2, 'MarkerSize', 8, ...
         'DisplayName', 'Pure FEM tot')
    %plot(margin_BEM_coeff, t_FEM_BEM, 'b*', 'MarkerSize', 12, 'LineWidth', 2, ...
    %     'DisplayName', 'FEM-BEM solve')
    plot(margin_BEM_coeff, t_FEM_BEM+TIME_mesh_bem_generation, 'r*', 'MarkerSize', 12, 'LineWidth', 2, ...
         'DisplayName', 'FEM-BEM tot')
    xlabel('Domain margin / Ly');
    ylabel('Wall-clock time (s)');
    title('Time comparison');
    legend; grid on

    salva_fig('time', folder_name, SALVATAGGIO);
    
    % --- Tabella riassuntiva ---
    fprintf('\n%-15s %-10s %-12s %-12s %-10s\n', ...
        'Margin/Ly', 'N nodes', 'err_L2 (%)', 'err_C_all', 'time (s)');
    fprintf('%-15s %-10s %-12s %-12s %-10s\n', ...
        '---------', '-------', '------', '------', '--------');
    for k = 1:n_test
        fprintf('%-15.0f %-10d %-12.4g %-12.4g %-10.3g\n', ...
            margini_label(k), N_nodes(k), err_V(k)*100, err_C_energy_all(k), t_fem_solve(k)+t_fem_mesh(k));
    end
    fprintf('%-15s %-10d %-12s %-12.4g %-10.3f\n', ...
        'FEM-BEM(ref)', length(x_r), '--', C_energy_all, t_FEM_BEM+TIME_mesh_bem_generation);


    if(SALVATAGGIO == 1)
        %Salva le variabili
        save(fullfile(folder_name, 'workspace_error.mat'));
    end
        
        
end