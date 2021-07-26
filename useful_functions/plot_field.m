function [fig] = plot_field(fig,field_data,plot_data,fonts_data)

name       = field_data.name;   % field name
y          = field_data.y;      % field data


flag = 1; % 1 if transient

% check if data refers to steady-state or transient problem
if ~isfield(plot_data,'tt')
    tt   = 1;
    flag = 0;
else
    tt     = plot_data.tt;
    dt     = plot_data.dt;
    time = dt*(tt-1);
end

limits = plot_data.limits;
plot_title  = plot_data.title;

font_title = fonts_data.font_title;
font_label = fonts_data.font_label;

if isfield(field_data,'reduced')
    
    field_vertices = field_data.reduced.vertices;
    field_elements = field_data.reduced.elements;
    field_indexes  = field_data.reduced.indexes;
else
    mesh       = field_data.mesh;
    field_vertices = mesh.vertices;
    field_elements = mesh.elements;
    field_indexes  = 1:length(mesh.vertices);
end


%Plot field
fig(length(fig)).Name = strcat(name,"_field_",num2str(tt));
pdeplot(field_vertices,field_elements(1:3,:),'xydata',y(field_indexes,tt),'xystyle','interp',...
        'zdata',y(field_indexes,tt)+2*abs(min(y(field_indexes,tt))) ,'zstyle','continuous',...
       'colorbar','on','mesh','off')

if flag
    title_temp = title(strcat(plot_title,' field  - Time = $',num2str(time,'%4.2f'),'$ [s]'),'interpreter','latex');
else
    title_temp = title(strcat(plot_title,' field '),'interpreter','latex'); 
end


set(title_temp,'FontSize',font_title);
label_temp = xlabel('$x_1$','interpreter','latex');
set(label_temp,'FontSize',font_label);
label_temp = ylabel('$x_2$','interpreter','latex');
set(label_temp,'FontSize',font_label);

% set domain limits
axis equal
if ~isfield(plot_data,'dom_limits')
    axis([-1.1 1.1 -1.1 1.1]);
else
    axis(plot_data.dom_limits);
end

colormap(jet);
caxis(limits)
view(0,90)
colorbar;
    
    
end

