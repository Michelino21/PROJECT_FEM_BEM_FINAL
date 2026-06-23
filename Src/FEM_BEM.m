%% FEM_BEM.m
    
function FEM_BEM(params)
    
    % Coupling FEM-BEM per condensatore a facce piane parallele 2D
    %addpath(genpath('.'));
    
    SALVATAGGIO      = params.SALVATAGGIO;
    ERROR_ANALYSIS   = params.ERROR_ANALYSIS;

    folder_name = '';
    
    %% SALVATAGGI
    if(SALVATAGGIO == 1)
    
        %nome cartella dinamico
        
        %delx
        % Moltiplica delx per 10000 e arrotonda per ottenere un intero di 4 cifre
        mesh_int = floor(params.delx * 1e4);
        
        % Controlla che il valore non superi 4 cifre (max 9999), altrimenti blocca
        if mesh_int > 9999
            error('delx = %.6f is too large: mesh index %d exceeds 4 digits.', params.delx, mesh_int);
        end
        
        % Controlla che il valore non sia troppo piccolo (minimo delx = 0.0001)
        if mesh_int < 1
            error('delx = %.6f is too small: minimum allowed value is 0.0001.', params.delx);
        end
        
        % Converte l'intero in stringa con zero padding a 4 cifre (es. 50 -> '0050')
        mesh_str = sprintf('%04d', mesh_int);

        %LY
        % Moltiplica delx per 10000 e arrotonda per ottenere un intero di 4 cifre
        Ly_int = floor(params.Ly * 1e4);
        
        % Controlla che il valore non superi 4 cifre (max 9999), altrimenti blocca
        if Ly_int > 9999
            error('Ly = %.6f is too large: mesh index %d exceeds 4 digits.', params.Ly, Ly_int);
        end
        
        % Controlla che il valore non sia troppo piccolo (minimo Ly = 0.0001)
        if Ly_int < 1
            error('Ly = %.6f is too small: minimum allowed value is 0.0001.', params.Ly);
        end
        
        % Converte l'intero in stringa con zero padding a 4 cifre (es. 50 -> '0050')
        Ly_str = sprintf('%04d', Ly_int);
        

        % LX
        % Moltiplica delx per 10000 e arrotonda per ottenere un intero di 4 cifre
        Lx_int = floor(params.Lx * 1e4);
        
        % Controlla che il valore non superi 4 cifre (max 9999), altrimenti blocca
        if Lx_int > 9999
            error('Lx = %.6f is too large: mesh index %d exceeds 4 digits.', params.Lx, Lx_int);
        end
        
        % Controlla che il valore non sia troppo piccolo (minimo Lx = 0.0001)
        if Lx_int < 1
            error('Lx = %.6f is too small: minimum allowed value is 0.0001.', params.Lx);
        end
        
        % Converte l'intero in stringa con zero padding a 4 cifre (es. 50 -> '0050')
        Lx_str = sprintf('%04d', Lx_int);
        % Costruisce il path della cartella con mesh e margine
        %folder_name = sprintf('./Documentation/runs/mesh_%s/margin_%d', mesh_str, params.margin_BEM_coeff);
        
        t = datetime('now','Format','yy_MM_dd__HH_mm_ss');

        folder_name = sprintf('./Documentation/runs/Lx_%s_Ly_%s/mesh_%s/margin_%d/%s', ...
            Lx_str, Ly_str, mesh_str, params.margin_BEM_coeff, string(t));

        folder_png  = fullfile(folder_name, 'png');
        folder_fig  = fullfile(folder_name, 'fig');
        
        % Crea la cartella se non esiste
        if ~exist(folder_name, 'dir')
            mkdir(folder_name);
            mkdir(folder_png);
            mkdir(folder_fig);
        end
        
        % Nome del file di log
        log_file = fullfile(folder_name, 'matlab_log.txt');
    
        % se esiste, cancellalo
        if exist(log_file, 'file')
            delete(log_file);
        end
        
        % Attiva il salvataggio della console
        diary(log_file)
        diary on
    
    end
    %% === PROBLEM PARAMETERS ===
    % --- Estrai parametri ---
    Lx               = params.Lx;
    Ly               = params.Ly;
    delx             = params.delx;
    dely             = params.dely;
    margin_BEM_coeff = params.margin_BEM_coeff;
    e0               = params.e0;
    er               = params.er;
    V1               = params.V1;
    V2               = params.V2;
    
    
    %parametri derivati
    margine = margin_BEM_coeff * Ly;     % = 0.03 m per lato
    X_mesh = Lx + 2*margine;  % = 0.16 m
    Y_mesh = Ly + 2*margine;  % = 0.07 m
    
    
    % ===== DISPLAY =====
    fprintf('\n--- GEOMETRIC PARAMETERS ---\n');
    fprintf('Lx               : %.4f m\n', Lx);
    fprintf('Ly               : %.4f m\n', Ly);
    fprintf('delx             : %.4f m\n', delx);
    fprintf('dely             : %.4f m\n', dely);
    
    fprintf('\n--- DOMAIN DIMENSIONS ---\n');
    fprintf('margin_BEM_coeff : %d\n',     margin_BEM_coeff);
    fprintf('margin           : %.4f m\n', margine);
    fprintf('X_mesh           : %.4f m\n', X_mesh);
    fprintf('Y_mesh           : %.4f m\n', Y_mesh);
    
    fprintf('\n--- PHYSICAL PARAMETERS ---\n');
    fprintf('e0               : %.4e F/m\n', e0);
    fprintf('er               : %.4f\n',     er);
    fprintf('V1               : %.4f V\n',   V1);
    fprintf('V2               : %.4f V\n',   V2);
    
    %% === MESH ===
    tic
    [conn, x, y, xmid, ymid, outb_node, ptop_node, pbottom_node, d_elm] = ...
        meshgen_cap(X_mesh, Y_mesh, Lx, Ly, delx, dely);    
    TIME_mesh_bem_generation = toc;
    
    
    M = size(conn, 1);
    N = length(x);
    
    eps_e = ones(M,1) * e0;
    eps_e(d_elm) = er * e0;
    
    fprintf('\n--- MESH ---\n');
    fprintf('Time of mesh generation: %.4f s\n', TIME_mesh_bem_generation);
    fprintf('Mesh Dimension: %d\n', N);
    fprintf('Connectivity Dimension: %d\n', M);
    
    
    %% === PARTIZIONE DEI NODI ===
    D_nodes = [ptop_node;    pbottom_node]; %Plates Nodes (Dirichlet Condition)
    B_nodes = outb_node;                    %Boundary Nodes
    F_nodes = setdiff(setdiff(1:N, D_nodes), B_nodes)'; %All the other nodes
    
    number_of_B_nodes = length(B_nodes);
    fprintf('\n--- NODES ---\n');
    fprintf('Nodes on the plates: %d\n', length(D_nodes));
    fprintf('Nodes on the boundary: %d\n', number_of_B_nodes);
    fprintf('Other Nodes: %d\n', length(F_nodes));
    
    %% === ASSEMBLAGGIO S (matrice FEM globale) ===
    tic
    S = sparse(N, N);
    area = zeros(M,1);
    
    for e = 1:M
        Ne = 3; %Number of edges of the element, 3 is triangular!
        [Ae, area(e)] = assemble_S(e, x, y, conn, eps_e(e), eps_e(e), 0);
        
        for i = 1:Ne
            ig = conn(e,i);
            for j = 1:Ne
                jg = conn(e,j);
                S(ig,jg) = S(ig,jg) + Ae(i,j);
            end
        end
    end
    TIME_building_S = toc;
    
    nnz_el_S = nnz(S);
    
    fprintf('\n--- ASSEMBLE S ---\n');
    fprintf('Time of assemble S: %.4f s\n', TIME_building_S);
    fprintf('Element in S (sparse): %d\n', nnz_el_S);
    fprintf('Sparsity Percentage: %.4g\n', nnz_el_S/(N*N)*100);
    
    
    
    %% === BLOCCHI DI S ===
    S_FF = S(F_nodes, F_nodes);
    S_FB = S(F_nodes, B_nodes);
    S_FD = S(F_nodes, D_nodes); %RHS
    S_BF = S(B_nodes, F_nodes);
    S_BB = S(B_nodes, B_nodes);
    S_BD = S(B_nodes, D_nodes); %RHS
    
    %% === RHS ===
    u_D = [V1 * ones(length(ptop_node),    1);
           V2 * ones(length(pbottom_node),  1)];
    
    %% === ASSEMBLAGGIO T_BB ===
    tic
    T_BB = assemble_T(B_nodes, x, y);
    T_BB = e0 * T_BB; % anche questa parte di equazione è moltiplicata per e0
    TIME_building_T = toc;
    
    nnz_el_T = nnz(T_BB);
    
    fprintf('\n--- ASSEMBLE T ---\n');
    fprintf('Time of assemble T: %.4f s\n', TIME_building_T);
    fprintf('Element in T (sparse): %d\n', nnz_el_T);
    fprintf('Sparsity Percentage: %.4g\n', nnz_el_T/(number_of_B_nodes*number_of_B_nodes)*100);
    
    %% === BEM ===
    % Riordina i nodi di Gamma come percorso chiuso antiorario
    B_ordered = reorder_boundary(B_nodes, x, y);
    
    % Estrai coordinate nel giusto ordine per il BEM
    xB = x(B_ordered);
    yB = y(B_ordered);
    
    tic
    [G_bem, H_bem] = assemble_G_H(1, xB, yB); %1 means linear element
    TIME_assemble_G_H = toc;
    
    % permutazione per ordinare G\H con i nodi nell'ordine FEM e non ANTIORARIO
    [~, perm] = ismember(B_nodes, B_ordered);
    
    % permutazione per ordinare i nodi dall'ordine FEM all'ANTIORARIO
    [~, inv_perm] = ismember(B_ordered, B_nodes);
    
    number_of_el_G = numel(G_bem);
    number_of_el_H = numel(H_bem);
    fprintf('\n--- ASSEMBLE G H ---\n');
    fprintf('Time of assemble G and H: %.4f s\n', TIME_assemble_G_H);
    
    fprintf('Element in G: %d\n', number_of_el_G);
    fprintf('Sparsity Percentage: %.4f\n', nnz(G_bem)/number_of_el_G*100);
    fprintf('cond(G_bem) = %.4g\n', cond(full(G_bem)));
    fprintf('Lower Eigenvalue of G %.4g\n',min(eig(G_bem)));
    fprintf('\n');
    fprintf('Element in H: %d\n', number_of_el_H);
    fprintf('Sparsity Percentage: %.4f\n', nnz(H_bem)/number_of_el_H*100);
    fprintf('cond(H_bem) = %.4g\n', cond(full(H_bem)));
    fprintf('Diagonal values of H: %s\n', num2str(uniquetol(diag(H_bem), 1e-6)', '%.4f  '));
    
    %% === SISTEMA ACCOPPIATO ===
    
    % Solve G^-1 * H
    tic
    BEM_matrix = G_bem \ H_bem;
    TIME_solve_BEM_MATRIX = toc;
    
    BEM_matrix_reordered = BEM_matrix(perm,perm);
    coupling_term = T_BB * BEM_matrix_reordered;
    
    % Coupled System
    K = [S_FF,                        S_FB; ...
         S_BF,   S_BB + coupling_term];
    
    nnz_el_K = nnz(K);
    number_of_el_K = numel(K);
    
    rhs = -[S_FD; S_BD] * u_D;
    
    nnz_el_RHS = nnz(rhs);
    number_of_el_RHS = length(rhs);
    
    fprintf('\n--- COUPLED SYSTEM ---\n');
    fprintf('Time of solve BEM Matrix: %.4f s\n', TIME_solve_BEM_MATRIX);
    fprintf('cond(BEM_matrix) = %.4g\n', cond(full(BEM_matrix_reordered)));
    fprintf('cond(T_BB * BEM_matrix) = %.4g\n', cond(full(coupling_term)));
    fprintf('\n');
    fprintf('Number of elements of K: [%d x %d]\n', size(K, 1), size(K, 2));
    fprintf('Element of matrix K %d\n', number_of_el_K)
    fprintf('Element in K (semi-sparse): %d\n', nnz_el_K);
    fprintf('Sparsity Percentage: %.4g\n', nnz_el_K/number_of_el_K*100);
    fprintf('\n');
    fprintf('Element RHS: %d\n', number_of_el_RHS); 
    fprintf('Sparsity Percentage: %.4g\n', nnz_el_RHS/number_of_el_RHS*100);
    
    figure;
    imagesc(full(BEM_matrix));
    colorbar;
    colormap(jet);
    title('Heatmap of BEM matrix (G\\H)');
    xlabel('Columns');
    ylabel('Rows');
    axis equal tight;
    salva_fig('heatmap_BEM_matrix', folder_name, SALVATAGGIO);
    
    figure;
    spy(K);
    title('Heatmap of system matrix (K)');
    xlabel('Columns');
    ylabel('Rows');
    salva_fig('heatmap_K_matrix', folder_name, SALVATAGGIO);
    
    %% === SOLUZIONE ===
    tic
    sol   = K \ rhs;
    TIME_solve_SYSTEM = toc;
    
    fprintf('\n--- SOLVE COUPLED SYSYEM ---\n');
    fprintf('Time of solve COUPLED SYSTEM: %.4f s\n', TIME_solve_SYSTEM);
    
    t_FEM_BEM = TIME_building_S + TIME_building_T + TIME_assemble_G_H + TIME_solve_BEM_MATRIX + TIME_solve_SYSTEM;
    fprintf('Time FEM-BEM method (except Mesh Generation): %.4f s\n', t_FEM_BEM);
    
    u_F   = sol(1:length(F_nodes));
    u_B   = sol(length(F_nodes)+1:end);
    
    % Ricostruzione del potenziale globale (Step 1)
    u = zeros(N, 1);
    u(F_nodes) = u_F;
    u(B_nodes) = u_B;
    u(D_nodes) = u_D;
    
    fprintf('\n\n\n');
    %% ====== POST-PROCESSING ===================
    %
    %============================================
    fprintf('\n==========================================\n');
    fprintf('            POST-PROCESSING               \n');
    fprintf('==========================================\n\n');
    % Step 1.a: Potenziale Interno
    % Step 1.b: Campo Elettrico Interno e Analisi elementi di campo elettrico problematici
    % Step 2: Derivata Normale su Gamma
    % Step 3: Potenziale Esterno
    % Step 4.a: Densita superficiale di carica
    % Step 4.b: Capacità Carica Piastre
    % Step 4.c: Capacità Energia Dominio FEM
    
    %% === POTENZIALE INTERNO === (Step 1.a)
    
    % --- Plot potenziale (mappa colori) ---
    [up, xp, yp] = create2darray(x, y, u, delx, dely);
    
    figure('Color', 'w');  % crea la figura con sfondo bianco
    clf;                   % pulisce la figura
    imagesc(xp, yp, up);
    colormap(jet); axis equal tight;
    set(gca,'Ydir','normal'); colorbar;
    xlabel('x (m)'); ylabel('y (m)');
    title('Potential (V)');
    set(gcf,'Color',[1 1 1]); hold on
    plot(x(ptop_node),    y(ptop_node),    'k', 'linewidth', 5)
    plot(x(pbottom_node), y(pbottom_node), 'k', 'linewidth', 5)
    
    salva_fig('potential_map', folder_name, SALVATAGGIO);
    
    % --- Plot potenziale (linee di livello) ---
    figure('Color', 'w');  % crea la figura con sfondo bianco
    clf;                    % pulisce la figura
    contour(xp, yp, up, 'LevelStep', 0.1)
    colormap(jet); axis equal tight; colorbar;
    xlabel('x (m)'); ylabel('y (m)');
    title('Potential (V) - contour');
    set(gcf,'Color',[1 1 1]); hold on
    plot(x(ptop_node),    y(ptop_node),    'k', 'linewidth', 5)
    plot(x(pbottom_node), y(pbottom_node), 'k', 'linewidth', 5)
    
    salva_fig('potential_contour', folder_name, SALVATAGGIO);
    
    %% === CAMPO ELETTRICO INTERNO === (Step 1.b)
    
    % --- Campo elettrico per elemento ---
    Ex = zeros(M,1); Ey = zeros(M,1);
    for e = 1:M
        [Ex(e), Ey(e)] = Efield(e, x, y, conn, u);
    end
    Emag = sqrt(Ex.^2 + Ey.^2);
    
    % --- Plot campo elettrico ---
    [Emagp, xp, yp] = create2darray(xmid, ymid, Emag, delx, dely);
    
    figure('Color', 'w');  % crea la figura con sfondo bianco
    clf;                    % pulisce la figura
    imagesc(xp, yp, Emagp);
    colormap(jet); axis equal tight;
    set(gca,'Ydir','normal'); colorbar;
    xlabel('x (m)'); ylabel('y (m)');
    title('Electric Field Intensity (V/m)');
    set(gcf,'Color',[1 1 1]); hold on
    quiver(xmid, ymid, Ex, Ey, 'color', 'k')
    plot(x(ptop_node),    y(ptop_node),    'k', 'linewidth', 5)
    plot(x(pbottom_node), y(pbottom_node), 'k', 'linewidth', 5)
    
    salva_fig('Efield_map', folder_name, SALVATAGGIO);
    
    
    % Visualizza distribuzione di E su d_elm
    figure
    histogram(Emag(d_elm), 50)
    xlabel('|E| (V/m)'); ylabel('Conteggio elementi');
    title('Electric Field distribution')
    
    salva_fig('Efield_histogram', folder_name, SALVATAGGIO);
    
    % Trova gli elementi problematici
    fprintf('------- E Field -------\n')
    E_teorico = (V1 - V2) / Ly;
    fattore_soglia = 5; % in percentuale
    soglia = (1 + fattore_soglia/100) * E_teorico;
    bad_elm  = d_elm(Emag(d_elm) > soglia);
    good_elm = d_elm(Emag(d_elm) <= soglia);
    
    
    fprintf('Theoretical E           = %.4f V/m\n', E_teorico);
    fprintf('Threshold Factor        = %.4f %%\n', fattore_soglia);
    fprintf('Threshold E             = %.4f V/m\n', soglia);
    fprintf('Elements with E > thr   : %d\n', length(bad_elm));
    fprintf('Elements with E <= thr  : %d\n', length(good_elm));
    
    % Visualizza dove sono i bad elements
    
    figure
    triplot(conn, x, y, 'Color', [0.8 0.8 0.8]); hold on
    if(~isempty(bad_elm))
        triplot(conn(bad_elm,:), x, y, 'r');
    end
    plot(x(ptop_node), y(ptop_node), 'k', 'linewidth', 3)
    plot(x(pbottom_node), y(pbottom_node), 'k', 'linewidth', 3)
    title('Elementi con E anomalo (rosso)')
    
    salva_fig('Efield_bad_elements', folder_name, SALVATAGGIO);
    
    
    %% === DERIVATA NORMALE SU GAMMA (per post-processing e soluzione esterna) === (Step 2)
    % Segno più perché q_FEM = q_BEM
    q_FEM = (BEM_matrix_reordered) * u_B;   % ∂u/∂n con normale uscente da Omega_int
    
    % Riordina q_FEM in senso antiorario.
    q_ordered = q_FEM(inv_perm);
    
    
    % grafico
    tol_idx = 1e-6;
    idx_bot   = find(abs(yB - min(yB)) < tol_idx);  % tutti i nodi con y = y_min
    idx_top   = find(abs(yB - max(yB)) < tol_idx);  % tutti i nodi con y = y_max
    idx_left  = find(abs(xB - min(xB)) < tol_idx);  % tutti i nodi con x = x_min
    idx_right = find(abs(xB - max(xB)) < tol_idx);  % tutti i nodi con x = x_max
    
    % --- Colori di sfondo per i 4 lati ---
    colors_patch = [1.0  0.85 0.85;   % Bottom - rosso chiaro
                    0.85 1.0  0.85;   % Right  - verde chiaro
                    0.85 0.90 1.0;    % Top    - blu chiaro
                    1.0  1.0  0.80];  % Left   - giallo chiaro
    labels    = {'Bottom', 'Right', 'Top', 'Left'};
    idx_sides = {idx_bot(2:end-1), idx_right(2:end-1), idx_top(2:end-1), idx_left(2:end-1)};
    
    figure;
    set(gcf, 'Color', 'w');
    ax1 = axes(); hold on; grid on;
    
    y_min = min(q_ordered)*1.3;
    y_max = max(q_ordered)*1.3;
    
    % Bande verticali per i 4 lati
    for s = 1:4
        idx_s   = idx_sides{s};
        x_start = idx_s(1);
        x_end   = idx_s(end);
        patch([x_start x_end x_end x_start], ...
              [y_min y_min y_max y_max], ...
              colors_patch(s,:), 'EdgeColor','none', 'FaceAlpha', 1);
        text(mean(idx_s), y_max*0.92, labels{s}, ...
            'HorizontalAlignment','center', 'FontSize', 9, 'FontWeight','bold');
    end
    
    plot(1:length(q_ordered), q_ordered, 'k.-', 'MarkerSize', 12);
    yline(0, 'k--', 'LineWidth', 1.2);
    
    ylim([y_min, y_max]);
    xlabel('Node index on \Gamma');
    ylabel('\partial u / \partial n');
    title('$\frac{\partial u}{\partial n}$ on $\Gamma$ $(E \cdot \hat{n},\ \hat{n}\ \mathrm{outward})$', ...
        'Interpreter', 'latex');
    subtitle('q > 0 : outgoing  | q < 0 : incoming', ...
        'Interpreter', 'tex');
    
    salva_fig('q_gamma', folder_name, SALVATAGGIO);
    
    % --- Flusso di E
    fluxE = sum(q_FEM);
    fprintf('\n');
    fprintf('Flux of E on Gamma: %.4g V/m\n', fluxE);
    
    %% === DENSITA' DI CARICA SULLE PIASTRE  === (Step 4)
    
    % ========================================================================
    % Recupera E su piastra alta e bassa
    % ========================================================================
    
    
    % Trova elementi che contengono almeno un nodo della piastra superiore
    elements_top = [];
    elements_top_out = [];
    for e = 1:M
        nodes_in_element = conn(e, :);
        if any(ismember(nodes_in_element, ptop_node))
            
            
            y_centroid = mean(y(nodes_in_element));
            y_plate_top = y(ptop_node(1));  % y della piastra superiore

            % Verifica che l'elemento sia nel gap (sotto la piastra)
            if y_centroid < y_plate_top  % Elemento sotto la piastra
                elements_top = [elements_top; e];
            end

            % Verifica che l'elemento sia nel gap (sopra la piastra)
            if y_centroid > y_plate_top  % Elemento sotto la piastra
                elements_top_out = [elements_top_out; e];
            end


        end
    end

    % Riordina da sinistra a destra (x crescente del centroide)
    [~, ord] = sort(xmid(elements_top));
    elements_top     = elements_top(ord);
    [~, ord] = sort(xmid(elements_top_out));
    elements_top_out = elements_top_out(ord);

    % Trova elementi adiacenti alla piastra inferiore
    elements_bottom = [];
    elements_bottom_out = [];
    for e = 1:M
        nodes_in_element = conn(e, :);
        if any(ismember(nodes_in_element, pbottom_node))
            % Verifica che l'elemento sia nel gap (sopra la piastra)
            y_centroid = mean(y(nodes_in_element));
            y_plate_bottom = y(pbottom_node(1));  % y della piastra inferiore
            if y_centroid > y_plate_bottom  % Elemento sopra la piastra
                elements_bottom = [elements_bottom; e];
            end

            % Verifica che l'elemento sia nel gap (sotto la piastra)
            if y_centroid < y_plate_bottom  % Elemento sotto la piastra
                elements_bottom_out = [elements_bottom_out; e];
            end
        end
    end

    % Riordina da sinistra a destra (x crescente del centroide)
    [~, ord] = sort(xmid(elements_bottom));
    elements_bottom    = elements_bottom(ord);
    [~, ord] = sort(xmid(elements_bottom_out));
    elements_bottom_out = elements_bottom_out(ord);
    
    % ========================================================================
    % Fai la proiezione sulla normale
    % ========================================================================
    
    % Normale alla piastra superiore (verso il basso, verso il gap)
    n_top = [0; -1];
    
    % Normale alla piastra inferiore (verso l'alto, verso il gap)
    n_bottom = [0; 1];
    
    % Proiezione del campo E sulla normale per piastra superiore
    En_top = zeros(length(elements_top), 1);
    for i = 1:length(elements_top)
        e = elements_top(i);
        E_vec = [Ex(e); Ey(e)];
        En_top(i) = dot(E_vec, n_top);
    end

    % Proiezione del campo E sulla normale per piastra superiore (elementi
    % esterni)
    En_top_out = zeros(length(elements_top_out), 1);
    for i = 1:length(elements_top_out)
        e = elements_top_out(i);
        E_vec_out = [Ex(e); Ey(e)];
        En_top_out(i) = dot(E_vec_out, -n_top);
    end
    
    % Proiezione del campo E sulla normale per piastra inferiore
    En_bottom = zeros(length(elements_bottom), 1);
    for i = 1:length(elements_bottom)
        e = elements_bottom(i);
        E_vec = [Ex(e); Ey(e)];
        En_bottom(i) = dot(E_vec, n_bottom);
    end

    En_bottom_out = zeros(length(elements_bottom_out), 1);
    for i = 1:length(elements_bottom_out)
        e = elements_bottom_out(i);
        E_vec_out = [Ex(e); Ey(e)];
        En_bottom_out(i) = dot(E_vec_out, -n_bottom);
    end
    
    % ========================================================================
    % Calcola la densità di carica per elemento
    % ========================================================================
    
    % Densità di carica per elementi adiacenti alla piastra superiore
    sigma_top = e0 * er * En_top;

    % Densità di carica per elementi adiacenti alla piastra superiore OUT
    sigma_top_out = e0 * er * En_top_out;
    
    % Densità di carica per elementi adiacenti alla piastra inferiore
    sigma_bottom = e0 * er * En_bottom;

    sigma_bottom_out = e0 * er * En_bottom_out;
    
    % ========================================================================
    % Calcola Q (carica totale per unità di lunghezza)
    % ========================================================================
    
    % Per ogni elemento, calcola l'area (o lunghezza del lato sulla piastra)
    % Assumendo elementi triangolari in 2D
    
    % Carica sulla piastra superiore
    Q_top = 0;
    for i = 1:length(elements_top)
        %e = elements_top(i);
        %nodes_in_element = conn(e, :);
        
        % Trova il lato dell'elemento che giace sulla piastra
        % (i nodi che appartengono a ptop_node)
        %nodes_on_plate = nodes_in_element(ismember(nodes_in_element, ptop_node));
        


        Q_top = Q_top + sigma_top(i) * delx/2;
        

        % if length(nodes_on_plate) >= 2
        % 
        %     length_segment = abs(x(nodes_on_plate(2)) - x(nodes_on_plate(1)));
        % 
        %     % Contributo alla carica totale
        %     Q_top = Q_top + sigma_top(i) * length_segment;
        % end

    end

    Q_top_out = 0;
    for i = 1:length(elements_top_out)
        %e = elements_top_out(i);
        %nodes_in_element = conn(e, :);
        
        % Trova il lato dell'elemento che giace sulla piastra
        % (i nodi che appartengono a ptop_node)
        %nodes_on_plate = nodes_in_element(ismember(nodes_in_element, ptop_node));
        
        Q_top_out = Q_top_out + sigma_top_out(i) * delx/2;
        % 
        % if length(nodes_on_plate) >= 2
        % 
        %     length_segment = abs(x(nodes_on_plate(2)) - x(nodes_on_plate(1)));
        % 
        %     % Contributo alla carica totale
        %     Q_top_out = Q_top_out + sigma_top_out(i) * length_segment;
        % end
    end
    
    % Carica sulla piastra inferiore
    Q_bottom = 0;
    for i = 1:length(elements_bottom)
        %e = elements_bottom(i);
        %nodes_in_element = conn(e, :);
        
        % Trova il lato dell'elemento che giace sulla piastra
        %nodes_on_plate = nodes_in_element(ismember(nodes_in_element, pbottom_node));
        
        %if length(nodes_on_plate) >= 2
            % Calcola lunghezza del lato sulla piastra
        %    x1 = x(nodes_on_plate(1));
        %    y1 = y(nodes_on_plate(1));
        %    x2 = x(nodes_on_plate(2));
        %    y2 = y(nodes_on_plate(2));
            
        %    length_segment = sqrt((x2-x1)^2 + (y2-y1)^2);
            
            % Contributo alla carica totale
        %    Q_bottom = Q_bottom + sigma_bottom(i) * length_segment;
        %end
        Q_bottom = Q_bottom + sigma_bottom(i) * delx/2;
    end
    
    Q_bottom_out = 0;
    for i = 1:length(elements_bottom_out)
        %e = elements_bottom_out(i);
        %nodes_in_element = conn(e, :);
        
        % Trova il lato dell'elemento che giace sulla piastra
        % (i nodi che appartengono a pbottom_node)
        %nodes_on_plate = nodes_in_element(ismember(nodes_in_element, pbottom_node));
        
        %if length(nodes_on_plate) >= 2
            % Calcola lunghezza del lato sulla piastra
            %x1 = x(nodes_on_plate(1));
            %y1 = y(nodes_on_plate(1));
            %x2 = x(nodes_on_plate(2));
            %y2 = y(nodes_on_plate(2));
            
            %length_segment = sqrt((x2-x1)^2 + (y2-y1)^2);
            
        %    length_segment = abs(x(nodes_on_plate(2)) - x(nodes_on_plate(1)));
            % Contributo alla carica totale
        %    Q_bottom_out = Q_bottom_out + sigma_bottom_out(i) * length_segment;
        %end
        Q_bottom_out = Q_bottom_out + sigma_bottom_out(i) * delx/2;
    end
    % ========================================================================
    % Risultati
    % ========================================================================
    fprintf('\n');
    fprintf('------- CHARGE on PLATES -------\n')
    fprintf('Total charge upper plate  (Q_top)   : %.6e C/m\n', Q_top);
    fprintf('Total charge upper plate  (Q_top_out)   : %.6e C/m\n', Q_top_out);
    fprintf('Total charge lower plate  (Q_bottom): %.6e C/m\n', Q_bottom);
    fprintf('Total charge lower plate  (Q_bottom_out)   : %.6e C/m\n', Q_bottom_out);
    fprintf('Sum (Q_tot)               : %.6e C/m\n', Q_top + Q_bottom + Q_top_out + Q_bottom_out);
    fprintf('Ratio Q_top/Q_bot         : %.6e\n', (Q_top+Q_top_out)/(Q_bottom+Q_bottom_out));
    
    % Densità di carica media
    sigma_top_mean = mean(sigma_top + sigma_top_out);
    sigma_bottom_mean = mean(sigma_bottom + sigma_bottom_out);
    
    fprintf('\n');
    fprintf('Average charge density upper plate: %.6e C/m^2\n', sigma_top_mean);
    fprintf('Average charge density lower plate: %.6e C/m^2\n', sigma_bottom_mean);
    
    % visualizza la distribuzione di carica
    figure;
    subplot(2,1,1);
    plot(sigma_top+sigma_top_out, 'o-');
    title('Charge Density - Upper Plate');
    xlabel('Element index');
    ylabel('\sigma [C/m^2]');
    grid on;
    
    subplot(2,1,2);
    plot(sigma_bottom+sigma_bottom_out, 'o-');
    title('Charge Density - Lower Plate');
    xlabel('Element index');
    ylabel('\sigma [C/m^2]');
    grid on;
    salva_fig('charge_density', folder_name, SALVATAGGIO);
    
    %%  === CAPACITA' Carica Piastre === (Step 4.b)
    
    % Differenza di potenziale
    Delta_V = V1 - V2;
    
    % Capacitanza usando la carica dalla piastra superiore
    C_from_top = abs(Q_top) / Delta_V;
    C_from_top_out = abs(Q_top_out) / Delta_V;
    
    % Capacitanza usando la carica dalla piastra inferiore
    C_from_bottom = abs(Q_bottom) / Delta_V;
    C_from_bottom_out = abs(Q_bottom_out) / Delta_V;
    
    % Capacitanza media (dovrebbero essere quasi uguali)
    C = (C_from_top + C_from_top_out + C_from_bottom + C_from_bottom_out) / 2;
    
    % ========================================================================
    % Confronto con soluzione analitica (condensatore piano ideale)
    % ========================================================================
    
    % Capacitanza teorica per condensatore piano infinito
    C_theoretical = e0 * er * Lx/Ly;
    
    % Errore percentuale
    error_percent = abs(C - C_theoretical) / C_theoretical * 100;

    fprintf('\n');
    fprintf('------- CAPACITANCE CHARGE METHOD -------\n');
    %fprintf('Capacitance (from upper plate)          : %.6e F/m\n', C_from_top);
    fprintf('Capacitance (from upper plate)          : %.6e F/m\n', C_from_top + C_from_top_out);
    %fprintf('Capacitance (from lower plate)          : %.6e F/m\n', C_from_bottom);
    fprintf('Capacitance (from lower plate)          : %.6e F/m\n', C_from_bottom + C_from_bottom_out);
    fprintf('Average capacitance                     : %.6e F/m\n', C);
    fprintf('Theoretical capacitance (infinite plate): %.6e F/m\n', C_theoretical);
    fprintf('Relative error                          : %.4f%%\n', error_percent);
    fprintf('\n');
    fprintf('Lx : %.6f m\n', Lx);
    fprintf('Ly = %.6f m\n', Ly);
    fprintf('Q_{top}    : %+.6e C/m\n', Q_top + Q_top_out);
    fprintf('Q_{bottom} : %+.6e C/m\n', Q_bottom + Q_bottom_out);
    fprintf('Q_{sum}    : %+.6e C/m\n', Q_top +  Q_top_out + Q_bottom + Q_bottom_out);
    fprintf('ΔV         : %.6f V\n', Delta_V);

    %%  === CAPACITA' Energia Dominio FEM === (Step 4.c)

    [C_energy_all,C_energy_inside,C_energy_clean] = compute_capacitance(V1,V2,e0,er,d_elm,conn,x,y,u,soglia,area,C_theoretical);

    
    %% ERROR
    if(ERROR_ANALYSIS == 1)
        margini_list = params.margini_list;
        error_FEM(Lx,Ly,delx,dely,e0,er,V1,V2,u,t_FEM_BEM,folder_name,SALVATAGGIO,margin_BEM_coeff,margini_list,x, y,outb_node, ptop_node, pbottom_node,TIME_mesh_bem_generation,soglia,C_theoretical,C_energy_all,C_energy_inside,C_energy_clean);
    end
    
    %% SALVATAGGI
    if(SALVATAGGIO == 1)
        %Salva le variabili
        save(fullfile(folder_name, 'workspace.mat'));
    
        % Fine del logging
        fprintf('\n\n\nEND OF LOG');
        diary off
    end
    
    
    
    fprintf('\nEND OF CODE');

end