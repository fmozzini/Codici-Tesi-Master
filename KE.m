function KEl=KEl(m,x,y,z,time)
% W: weigth
% Partendo da acc calcolo velocità (integro)
% trapz performs numerical integration via the trapezoidal method. 
% This method approximates the integration over an interval by breaking 
% the area down into trapezoids with more easily computable areas.
vel_x = cumtrapz(time,x);
vel_y = cumtrapz(time,y);
vel_z = cumtrapz(time,z);
vel = sqrt(vel_x.^2 + vel_y.^2 + vel_z.^2);
KEl = 0.5*m*vel.^2;
end 