%this program is calculate the transport curve
%how long will transport take
transp_duration = 6;
%how far are we transporting the atoms
transp_distance = 200;

starting_distance = 200;
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
%lens controller connection
%comm settings
port1 = 'COM6';
port2 = 'COM3';
ser1 = establish_connection(port1); %top lens
% ser2 = establish_connection(port2);

resetLens(ser1);
setPowerMode(ser1);
% 
tic
while toc < transp_duration
    x = position(toc);
    set = power(x)*1000;
    fastPower(ser1,set)
end
toc

% fastPower(ser1,0);
fclose(ser1);


function clearBytes(s)
disp('bytes available')
s.BytesAvailable
fread(s)
end

%put the optotune into power set mode
function setPowerMode(s)
clearBytes(s);
str = double('MwCA');
final = appendCRC16(str);
fwrite(s,final);
pause(0.2)
disp('set power mode')
s.BytesAvailable;
A = fread(s,s.BytesAvailable);
low = bin2dec([dec2bin(A(5),8) dec2bin(A(6),8)])
high = bin2dec([dec2bin(A(7),8) dec2bin(A(8),8)]);
high = typecast(uint16(high),'int16')
end

function resetLens(s)
clearBytes(s);
fprintf(s,'Start');
end

%setting focal power without any byteflushes or delays for speed
function fastPower(s,power)
start = double('PwDA');
if power > 6
    power = 6
end
if power < -6
    power = -6
end
set = power*200;
set = dec2bin(typecast(int16(set),'uint16'),16);
set1 = bin2dec(set(1:8));
set2 = bin2dec(set(9:16));
dummy1 = double('g');
dummy2 = double('g');
appended = appendCRC16([start set1 set2 dummy1 dummy2]);
fwrite(s,appended);
end

%specify the desired focal power (diopters) 
%and then set the optotune to it
function setPower(s,power)
disp('clearing bytes')
clearBytes(s);
start = double('PwDA');
if power > 10
    power = 10
end
if power < -10
    power = -10
end
set = power*200;
set = dec2bin(typecast(int16(set),'uint16'),16);
set1 = bin2dec(set(1:8));
set2 = bin2dec(set(9:16));
dummy1 = double('g');
dummy2 = double('g');
appended = appendCRC16([start set1 set2 dummy1 dummy2]);
fwrite(s,appended);
end

function myans=appendCRC16(byte_array)

[low,high]=CRC16(byte_array);
myans=[byte_array, low, high];

end

function [low_byte, hi_byte]=CRC16(decinput)
poly16='8005';
poly16=[1 convert2arr(dec2bin(hex2dec(poly16),16))];

poly=poly16;
register=zeros(1,length(poly)-1);

for kk=1:length(decinput)
    bin=convert2arr(dec2bin(decinput(kk),8));
    bin=flip(bin);
    message=bin;
    if kk==length(decinput)
       message=[message zeros(1,length(poly)-1)]; 
    end
    
    while ~isempty(message)
        pop=register(1);
        register=[register(2:end) message(1)];
        message=message(2:end);
        if pop==1
           register=register+poly(2:end);
           for jj=1:length(register)
              register(jj)= mod(register(jj),2);
           end
        end
    end
end

register=flip(register);
str=convert2str(register);
byte=bin2dec(str);
hi_byte=bin2dec(str(1:8));
low_byte=bin2dec(str(9:16));

function mystr=convert2str(arr)
mystr='';
for ii=1:length(arr)
   mystr=[mystr num2str(arr(ii))];
end
end

function arr=convert2arr(str)
arr=[];

for ii=1:length(str)
   arr=[arr str2num(str(ii))]; 
end
end

end

function [s] = establish_connection(PORT)
baud = 115200;
parity = 'none';
stop = 1;
data = 8;
terminator = 'CR/LF';
flowcontrol = 'none';
Timeout = 2;

s = serial(PORT);
fopen(s);
disp("successful connection estd.")
%set comm parameters
s.BaudRate = baud;
s.DataBits = data;
s.Parity = parity;
s.StopBits = stop;
s.Terminator = terminator;
s.flowcontrol = flowcontrol;
s.Timeout = Timeout;
end