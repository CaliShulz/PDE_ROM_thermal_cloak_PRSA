function [] = get_grey_scrofa(obs_vert)

%Generate shapes to fill domain not controlled



x_c_obs = obs_vert(1,:);
y_c_obs = obs_vert(2,:);

x_domain = [-1 -1 1  1 ] ;
y_domain = [-1  1 1 -1 ] ;

x = [ x_domain   nan x_c_obs ]; 
y = [ y_domain   nan y_c_obs ]; 

pgon_int = polyshape(x,y,'Simplify',true);
plot(pgon_int,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',1)



%a1 = alphaShape( pgon_vert(:,1), pgon_vert(:,2), -100*ones(size(pgon_vert,1), 1));

%plot(pgon_int,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',1)
%plot3(a1.Points(:,1),a1.Points(:,2),a1.Points(:,3),'b')


end



