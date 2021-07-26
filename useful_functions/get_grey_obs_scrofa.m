function [] = get_grey_obs_scrofa(outer_vert,obs_vert)

%Generate shapes to fill domain not controlled




x = [ outer_vert(1,:) nan obs_vert(1,:)]; 
y = [ outer_vert(2,:) nan obs_vert(2,:) ]; 
pgon_int = polyshape(x,y,'Simplify',true);


plot(pgon_int,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',1)



end



