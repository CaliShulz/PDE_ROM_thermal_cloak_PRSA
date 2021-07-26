function [] = get_grey_control(r_control_ext,r_control_in,r_obs)
%Generate shapes to fill domain not controlled


theta     = linspace(0,2*pi,100);
r_obs     = 0.2;

x_c_ext = r_control_ext * cos(theta);
y_c_ext = r_control_ext * sin(theta);

x_c_int = r_control_in * cos(theta);
y_c_int = r_control_in * sin(theta);

x_c_obs = r_obs * cos(theta);
y_c_obs = r_obs * sin(theta);



x_domain = [-1 -1 1  1] ;
y_domain = [-1  1 1 -1] ;


x = [ x_domain nan x_c_ext ]; 
y = [ y_domain nan y_c_ext ]; 
pgon_ext = polyshape(x,y,'Simplify',true);


x = [ x_c_int nan x_c_obs ]; 
y = [ y_c_int nan y_c_obs ]; 

pgon_int = polyshape(x,y,'Simplify',true);

plot(pgon_ext,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',1)
hold on
plot(pgon_int,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',1)



end

