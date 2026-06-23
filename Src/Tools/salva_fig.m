function salva_fig(nome, folder_name, SALVATAGGIO)
    if SALVATAGGIO == 1

        % FONT SIZE
        ax = gca;
        ax.FontSize = 20;

        saveas(gcf, fullfile(folder_name, 'fig', [nome '.fig']));
        exportgraphics(gcf, fullfile(folder_name, 'png', [nome '.png']), 'Resolution', 300);
    end
end