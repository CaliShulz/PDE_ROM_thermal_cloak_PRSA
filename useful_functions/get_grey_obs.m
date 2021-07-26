function [] = get_grey_obs(r_control_ext,r_obs)

%Generate shapes to fill domain not controlled


theta     = linspace(0,2*pi,100);

x_c_ext = r_control_ext * cos(theta);
y_c_ext = r_control_ext * sin(theta);


x_c_obs = r_obs * cos(theta);
y_c_obs = r_obs * sin(theta);


x = [ x_c_ext nan x_c_obs ]; 
y = [ y_c_ext nan y_c_obs ]; 
pgon_int = polyshape(x,y,'Simplify',true);


plot(pgon_int,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',1)



end



