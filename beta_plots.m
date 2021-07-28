
function [] = beta_plots()

mesh_names = dir('*.mat') ;
sslash = path_setup() ; % setup path 

eta_mesh = [];

fig = gobjects(0);
set(0,'DefaultFigureVisible','on');

fig(length(fig)+1) = figure;
fig(length(fig)).Name = "beta_comparison_ss";

for tt = 1:length(mesh_names)

    mesh_name = erase(mesh_names(tt).name,'.mat');
    
    [eta_plot_tt] = beta_comparison(mesh_name,[3.5 15000     0     1E-10]);
    
    semilogx(eta_plot_tt(2,:),eta_plot_tt(1,:),'o-','Linewidth',1.5);
    legend_name = erase(strrep(erase(mesh_name,'mesh_'),'_',' '),'cloak')
    Legends{tt} = legend_name;
    hold on;
    
end


font_label = 18;
font_title = 19;
font_legend = 10;


grid('minor');
set(gca, 'XDir','reverse');
ylim([0.965 1.015])
label_temp = xlabel('$\beta$','interpreter','latex');
set(label_temp,'FontSize',font_label);
label_temp = ylabel('$ \eta  $','interpreter','latex');
set(label_temp,'FontSize',font_label);
title_temp = title('Control weighting selection','interpreter','latex');
set(title_temp,'FontSize',font_title);
legend_temp = legend(Legends,'interpreter','latex','Location','northwest');
set(legend_temp,'FontSize',font_legend);
grid on
axis square

mkdir('archive_sim');
name_tmp_jj = strcat('archive_sim',sslash,fig(1).Name,'.png');
exportgraphics(fig(1),name_tmp_jj);

end
