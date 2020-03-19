%the purpose of this matlab script is to calculate numerically the required
%foci of the two tunable lenses in the telescope for a given beam waist and
%focus away from the final lens. Although I can find the focus and the
%waist in terms of the foci, I am not sure if it is possible to find an
%analytic expression going the other way. At least it is much more
%difficult. 
%this program is calculate the transport curve
%how long will transport take
transp_duration = 2; %seconds
%how far are we transporting the atoms
transp_distance = 300; %mm

starting_distance = 200; %mm
%using simple integration, find constant
%in kinematic equations to get above constants
c = 6*transp_distance./(transp_duration.^3);
%calculate velocity (parabolic)
velocity = @(t) c.*t.*(transp_duration-t);
%find position by integrating velocity
position = @(t) integral(velocity,0,t)+starting_distance;
%plot position v, and a
subplot(2,2,2);
time = linspace(0,transp_duration,1000);
v = arrayfun(velocity,time);
plot(time,v);
xlabel('time (s)')
ylabel('velocity (mm/s)')
title("Velocity Profile")
grid on;
subplot(2,2,1);
x = arrayfun(position,time);
plot(time,x);
xlabel('time (s)')
ylabel('focal position after last lens, (mm)')
title("Transport Profile")
%x = arrayfun(position,time)
grid on;
subplot(2,2,3);
step = transp_duration/100;
a = diff(v)/step;
a = a/1000;
a = a/9.8;
plot(time(1:length(a)),a)
xlabel('time (s)')
ylabel("acceleration (fraction of g (9.8m/s^2)")
title("acceleration profile")
grid on;
%now we need to input the optical parameters of the system
%focal length of the static system
f = 250;
%distance between the flexible lenses
% d = 100;
d = 90; %as of 7/26/2019
%image after the first lens (if collimated input beam, then this is the
%focus of the first lens. Lets assume this for now
i = 1000;

power = @(x) ((x-f)*(d-i)-f^2)/((f^2)*(i-d));
focal = @(x) ((f^2)*(i-d))/((x-f)*(d-i)-f^2);
%x = arrayfun(@(t) integral(@(a) velocity(a)
powers = arrayfun(power,x);
focals = arrayfun(focal,x);
powers = powers*1000;
subplot(2,2,4);
plot(time,powers);
title("optical power")
ylabel('Diopters')
xlabel('Time (s)')
grid on;