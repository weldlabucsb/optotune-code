%The purpose of this program is to be a PID control feedback system on the
%focal position of the lens
%Author:Max Prichard
%Modified: 6/26
%All units of distance are in mm
%All units of time are in s
%Note; it is possible to characterize the system using only the fundamnetal
%optical measurements of the system, for instance the focal lengths of the
%lenses and the distances between the optical elements, but I imagine that
%the better way to do this would be instead to measure the beam waist and
%then use that to determine some of the parameters of the system

%parameters:
f1 = 100;
f3 = 250;
o = 100;

%transport curve:
transp_duration = 0.5;
function v = velocity(t)
v = 1000*t*(transp_duration-t);

function x = position(t)
x = integral(velocity(x)




