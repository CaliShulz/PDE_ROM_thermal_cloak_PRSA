function [] = get_grey_plate(r_obs)

%Generate shapes to fill domain not controlled


theta     = linspace(0,2*pi,100);



x_c_obs = r_obs * cos(theta);
y_c_obs = r_obs * sin(theta);

x_domain = [-1 -1 1  1] ;
y_domain = [-1  1 1 -1] ;

x = [ x_domain nan x_c_obs ]; 
y = [ y_domain nan y_c_obs ]; 
pgon_int = polyshape(x,y,'Simplify',true);


plot(pgon_int,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',1)



end



