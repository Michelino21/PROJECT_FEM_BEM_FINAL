 %% === DENSITA' DI CARICA SULLE PIASTRE  === (Step 4)
 
 % Dopo aver usato il codice, se si ricarica l'ambiente matlab, si può
 % usare questo per calcolare la capacità con il metodo della carica,
 % sull'ultimo dominio disponibile, quindi il fem di dimensione massima (in
 % genere)


    % ========================================================================
    % Recupera E su piastra alta e bassa
    % ========================================================================
    
    
    M = size(conn_k, 1);
    conn = conn_k;
    y = y_k;
    x = x_k;
    u = u_fem_k;
    ptop_node = ptop_k;

    Ex = zeros(M,1); Ey = zeros(M,1);
    for e = 1:M
        [Ex(e), Ey(e)] = Efield(e, x, y, conn, u);
    end

    % Trova elementi che contengono almeno un nodo della piastra superiore
    elements_top = [];
    for e = 1:M
        nodes_in_element = conn_k(e, :);
        if any(ismember(nodes_in_element, ptop_k))
            % Verifica che l'elemento sia nel gap (sotto la piastra)
            y_centroid = mean(y(nodes_in_element));
            y_plate_top = y(ptop_k(1));  % y della piastra superiore
            if y_centroid < y_plate_top  % Elemento sotto la piastra
                elements_top = [elements_top; e];
            end
        end
    end
    
  
    
    % ========================================================================
    % Fai la proiezione sulla normale
    % ========================================================================
    
    % Normale alla piastra superiore (verso il basso, verso il gap)
    n_top = [0; -1];
    
    % Proiezione del campo E sulla normale per piastra superiore
    En_top = zeros(length(elements_top), 1);
    for i = 1:length(elements_top)
        e = elements_top(i);
        E_vec = [Ex(e); Ey(e)];
        En_top(i) = dot(E_vec, n_top);
    end
    

    % ========================================================================
    % Calcola la densità di carica per elemento
    % ========================================================================
    
    % Densità di carica per elementi adiacenti alla piastra superiore
    sigma_top = e0 * er * En_top;
    
    % ========================================================================
    % Calcola Q (carica totale per unità di lunghezza)
    % ========================================================================
    
    % Per ogni elemento, calcola l'area (o lunghezza del lato sulla piastra)
    % Assumendo elementi triangolari in 2D
    
    % Carica sulla piastra superiore
    Q_top = 0;
    for i = 1:length(elements_top)
        e = elements_top(i);
        nodes_in_element = conn_k(e, :);
        
        % Trova il lato dell'elemento che giace sulla piastra
        % (i nodi che appartengono a ptop_node)
        nodes_on_plate = nodes_in_element(ismember(nodes_in_element, ptop_k));
        
        if length(nodes_on_plate) >= 2
            % Calcola lunghezza del lato sulla piastra
            %x1 = x(nodes_on_plate(1));
            %y1 = y(nodes_on_plate(1));
            %x2 = x(nodes_on_plate(2));
            %y2 = y(nodes_on_plate(2));
            
            %length_segment = sqrt((x2-x1)^2 + (y2-y1)^2);
            
            length_segment = abs(x(nodes_on_plate(2)) - x(nodes_on_plate(1)));
            % Contributo alla carica totale
            Q_top = Q_top + sigma_top(i) * length_segment;
        end
    end
 
    
    % ========================================================================
    % Risultati
    % ========================================================================
    fprintf('\n');
    fprintf('------- CHARGE on PLATES -------\n')
    fprintf('Total charge upper plate  (Q_top)   : %.6e C/m\n', Q_top);
   
     %%  === CAPACITA' Carica Piastre === (Step 4.b)
    
     V1 = 0.05;
     V2 = 0.05;

    % Differenza di potenziale
    Delta_V = V1 - V2;
    
    % Capacitanza usando la carica dalla piastra superiore
    C_from_top = abs(Q_top) / Delta_V;