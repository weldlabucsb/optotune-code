(* ::Package:: *)

function OptotuneComm
delete(instrfindall);
imaqreset
daqreset

s_A=establishConnection(4,'A');
s_B=establishConnection(3,'B');
AO=generateAO;

handshake(s_A);
getID(s_A);
getFirmware(s_A);
disp(' ');

handshake(s_B);
getID(s_B);
getFirmware(s_B);
disp(' ');

setAnalog(AO,0,0);

getTemp(s_A)
getTemp(s_B)

try 
%    0. DATA LOADING: Here we load anything we need from previous runs that
%     we will use in lieu of taking new data.
%     
%   I. TEMPERATURE COMPENSATION AND VERIFICATION: PREPARING FOR TEMPERATURE
%   COMPENSATION
%     
%     1. Steady state temperature data collection. This will be used to
%       determine the current needed to reach desired temperature for
%       driftCalibration function.
%
%     atempdata=ssTemp(s_A,s_B,140,4000,10,50);
%     btempdata=ssTemp(s_B,s_A,140,4000,10,50);
%     
%     2. Collect drift data for each lens individually (temp/curr space)
%
%     adrift=driftCalibration(AO,s_A,s_B,atempdata,25,30,140,4000,21,50,0.1);
%     bdrift=driftCalibration(AO,s_B,s_A,btempdata,25,30,140,4000,21,50,0.1);
%
%     3. Verify that the current interpolation method (compCurrInterp) can
%       remove temperature dependence from optical behavior.
%
%     aTestMatrix=driftCalTest(AO,s_A,s_B,atempdata,adrift,25,30,140,4000,21,50,0.1);
%     bTestMatrix=driftCalTest(AO,s_B,s_A,btempdata,bdrift,25,30,140,4000,21,50,0.1);
%     
%   II. ACQUIRE (UN)COMPENSATED CALIBRATION DATA IN CURRENT-CURRENT SPACE
%     
%     1. Randomly collect spot size data in the I_A-I_B current space. This
%       data is used to see how the minimum line (focus line) changes as we
%       move through current space. This can be done using raw currents
%       (uncompensated) or temperature-corrected currents (compensated).
%
%     keyboard
%     uncompcal=getCalAnalog(AO,s_A,s_B,20,32,140,3500,50);
%     keyboard
%     compcal=getCompCalAnalog(AO,s_A,s_B,20,32,140,4000,50,adrift,bdrift);
%     uncompcal=getCalRand(s_A,s_B,20,30,140,4000,10);
%     compcal=getCompCalRand(s_A,s_B,20,30,140,4000,50,adrift,bdrift);
%     
%   III. CALIBRATION FIT FUNCTIONS: CURVE FITTING TO FOLLOW POINTS OF FOCUS
%
%     1. Here we look at the (un)compensated calibration data sets
%       collected in II. We follow the minimum line, generating a fit
%       function I_A(I_B) that follows the minima; given a requested I_B, the
%       fit function chooses I_A such that the resulting spot size is a
%       minimum.
%
%     [uncompfocusobject,~]=focusFit(uncompcal,140,4000);
%     [compfocusobject,~]=focusFit(compcal,140,4000);
%     
%   IV. CALIBRATION TESTING: UNCOMPENSATED AND COMPENSATED
%
%     1. Simple testing of (un)compensated calibration. Here we run through
%       the motions of following the minimum line without taking any data.
%       These functions are good for demonstrations.
% 
%     uncompCalTest(s_A,s_B,uncompfocusobject,imin,imax,steps,delay)
%     compCalTest: Needs to be remade following the transition from
%       compCurrent to compCurrInterp
%
%     2. Same as before but now we're taking data to see how well we're
%       following the intended focsued spot sizes.
%
%     uncompCalTestDataRand(s_A,s_B,uncompcal,uncompfocusobject,140,4000,30,0.25);
%     compCalTestDataRand(s_A,s_B,compcal,compfocusobject,adrift,bdrift,140,4000,30,0.25);

keyboard


resetCurr(s_A,s_B);
fclose(s_A);
fclose(s_B);

disp('Closing all connections.');

delete(instrfindall)
imaqreset
daqreset
    
catch
    disp('Exception found!');
    resetCurr(s_A,s_B);
    fclose(s_A);
    fclose(s_B);
    delete(instrfindall)
    imaqreset
end
end

% GENERAL SERIAL/ANALOG COMMANDS
function clearBytes(s)
% This function clears any serial communication bytes that are sitting in
% channel s
while s.BytesAvailable>0
    fread(s,1);
end
end

function s=establishConnection(portno,portname)
% This function establishes a serial connection between MATLAB and the
% desired com port
port=['COM' num2str(portno)];
baud=115200;
databits=8;
stopbits=1;
parity='none';
flowcontrol='none';

display(sprintf('Opening % s...',port));
s = serial(port,'BaudRate',baud,'Parity',parity,'StopBits',stopbits,'DataBits'
	,databits,'FlowControl',flowcontrol,'Name',portname);

set(s,'terminator','CR');
fopen(s);

disp(s);
disp('Serial connection established.');
disp(' ');

end

function AO=generateAO
% This function uses the MCC USB-1408FS-Plus and creates an analog output
% object (AO). It then adds the two analog output channels (0 and 1) and
% renames them a_A and a_B (analog output for driver A and B, respectively).
%The second parameter on analogoutput() might need to be changed depending
% on whether or not there are other MCC devices on the computer.

AO=analogoutput('mcc',1);
addchannel(AO,0:1,{'a_A','a_B'});
AO.SampleRate=25000;
disp('USB-1408FS-Plus: AO device ready.')
disp(' ');
end

% GENERAL READ/WRITE COMMANDS: Talk to the lens drivers and DAQ.
function handshake(s)
% Is the Lens Driver available?
clearBytes(s)
fprintf(s,'% s \r','Start');
reply=fread(s,7);
disp(['Lens ' get (s,'Name') ' ' char(transpose(reply(1:5)))]);
end

function getID(s)
% This gets the ID of the Lens Driver connected to port s. This also
% verifies that the driver is available.
clearBytes(s)
mystr='IRaaaaaaaa';%5993
bin=uint8(mystr);
high=uint8(hex2dec('59'));
low=uint8(hex2dec('93'));

msg=[bin low high];

fwrite(s,msg);
reply=fread(s,14);

disp(['Board ID ' get (s,'Name') ': ' char(transpose(reply(3:10)))]);
end

function getFirmware(s)
% Is the driver up to date?

clearBytes(s)
% mystr='Va';%48 FE
% bin=uint8(mystr);
% high=uint8(hex2dec('48'));
% low=uint8(hex2dec('FE'));
% msg=[bin low high];

fwrite(s,[86   128    62]);

reply=fread(s,11);
reply=reply(2:7);
a=num2str(reply(1));
b=num2str(reply(2));
c=num2str(bin2dec([dec2bin(reply(3)) dec2bin(reply(4))]));
d=num2str(bin2dec([dec2bin(reply(5)) dec2bin(reply(6))]));

disp(['Firmware Version ' get (s,'Name') ': ' a '.' b '.' c '.00' d]);
end

function getStatus(s)
clearBytes(s)

fwrite(s,[83   114   188   213])
reply=fread(s,9);

% sb4=dec2bin(reply(2),8);
% sb3=dec2bin(reply(3),8);
sb2=dec2bin(reply(4),8);
sb1=dec2bin(reply(5),8);

msg1={'- Temperature out of range specified by user.'
    '- Focal power out of guaranteed range (defined by user set temperature range.'
    '- Temperature is outside product specifications.'
    '- Focal power inversion.'
    '- Cannot reach lens focal power.'
    '- No temperature limits received for controlled mode.'
    '- No or faulty EEPROM.'
    '- Not all hardware available.'};
msg2='- The connected lens is not compatible with the firmware on the Lens Driver.';

disp(['Status ' get (s,'Name') ':'])

for ii=1:8
    if sb1(ii)==49
        disp(char(msg1(ii,:)));
    end
end

if sb2(8)==49
    disp(msg2)
end
end

function temp=getTemp(s)
% What is the temperature of the lens connected to the lens driver on port s?
clearBytes(s)
fwrite(s,[84    67    65   176   208]);
reply=fread(s,9);
temp=bin2dec ([dec2bin(reply(4),8) dec2bin(reply(5),8)])/16;
% disp(['Temperature ' get (s,'Name') ': ' num2str (temp) ' degC']);
% disp(' ');
end

function alog=getAnalog(s) % #ok<*DEFNU>
% This checks the ADC on the lens driver connected to s. The result is given
% as a 10 bit number 4*(0-1023), mapped from 0-5 V.
clearBytes(s)
% mystr=double('GAA');
% msg=appendCRC16(mystr);

fwrite(s,[71 65 65 64 117]);
reply=fread(s,7);
alog=4*bin2dec([dec2bin(reply(4),8),dec2bin(reply(5),8)]);
end

function volt=getVolt(s)
% Same as getAnalog(s) but now we're converting the measurement to volts.
clearBytes(s)
% mystr=double('GAA');
% msg=appendCRC16(mystr);

fwrite(s,[71 65 65 64 117]);
reply=fread(s,7);
out=bin2dec([dec2bin(reply(4),8),dec2bin(reply(5),8)]);
volt=5*1.0170*(out/1024);
end

function setAnalog(AO,atarget_new,btarget_new)
% This immediately pushes the 12 bit values atarget and btarget to
% channels a_A and a_B on the analog output device AO. Cal is:
% a_A: 5V in -> 991/1023=4.8436 V out
% a_B: 5V in -> 993/1023=4.8534 V out

if atarget_new > 4095
    atarget_new=4095;
end
if atarget_new < 0
    atarget_new=0;
end
if btarget_new > 4095
    btarget_new=4095;
end
if btarget_new < 0
    btarget_new=0;
end

putsample(AO,[5*atarget_new/4095 5*btarget_new/4095])
wait(AO,5)
end
function rampAnalog(AO,atarget_new,btarget_new)

if atarget_new > 4095
    atarget_new=4095;
end
if atarget_new < 0
    atarget_new=0;
end
if btarget_new > 4095
    btarget_new=4095;
end
if btarget_new < 0
    btarget_new=0;
end

persistent atarget_old
persistent btarget_old

if isempty(atarget_old)
    atarget_old=0;
    btarget_old=0;
end

asteps=5*transpose (linspace(atarget_old,atarget_new,100))/4095;
bsteps=5*transpose (linspace(btarget_old,btarget_new,100))/4095;
atarget_old=atarget_new;
btarget_old=btarget_new;

putdata(AO,[asteps bsteps]);
start(AO)
wait(AO,5)
end

% MODE COMMANDS: Change lens operating mode.
function getMode(s)
clearBytes(s)
fwrite(s,[77    77    65   101   119])
reply=fread(s,8);

modeArray={'DC'
    'Sinusoidal'
    'Triangular'
    'Rectangular'
    'Focal Power'
    'Analog'
    'Position Controlled'};

for ii=1:7
    if reply(4)==ii
        disp(['Mode ' get (s,'Name') ': ' char(modeArray(ii,:))]);
    end
end

disp(' ')
end

function modeDC(s)
% This function switches the lens to DC mode from any other mode. It seems
% to work as modeSerial is supposed to work when switching from Analog mode.

% mystr=double('MwDA');
% msg=appendCRC16(mystr);
fwrite(s,[77   119    68    65    84    70]);
clearBytes(s)
end

function modeAnalog(s)
% This switches the lens to Analog mode where the current is controlled by
% the 10 bit ADC on the lens driver.
fwrite(s,[77   119    65    65    87    22]);
clearBytes(s)
end

function modeSin(s)
% mystr=double('MwSA');
% msg=appendCRC16(mystr);
fwrite(s,[77   119    83    65    91   182]);
clearBytes(s);
end

function modeRec(s)
% mystr=double('MwQA');
% msg=appendCRC16(mystr);
fwrite(s,[77   119    81    65    90   214]);
clearBytes(s);
end

function modeTri(s)
% mystr=double('MwTA');
% msg=appendCRC16(mystr);
fwrite(s,[77   119    84    65    89   134]);
clearBytes(s);
end

function modeSerial(s)
% In theory this brings us back to serial control of the currents in the
% lenses, or at a minimum disables analog control. This might be useless as
% modeDC does a better job of serving this purpose.
fwrite(s,[76    83   116   253]);
clearBytes(s)
end

% PROPERTY READ/CHANGE COMMANDS: Manage operating mode parameters.
function high=getHigh(s)
% This gets the upper bound for current swings in the three periodic lens modes.
clearBytes(s);
% mystr=double('PrUAaaaa');
% msg=appendCRC16(mystr);
fwrite(s,[80   114    85    65    97    97    97    97   193   107]);
reply=fread(s,9);
high=bin2dec([dec2bin(reply(4),8) dec2bin(reply(5),8)]);
end

function low=getLow(s)
% This gets the lower bound for current swings in the three periodic lens modes.
clearBytes(s);
% mystr=double('PrLAaaaa');
% msg=appendCRC16(mystr);
fwrite(s,[80   114    76    65    97    97    97    97   195    98]);
reply=fread(s,9);
low=bin2dec([dec2bin(reply(4),8) dec2bin(reply(5),8)]);
end

function freq=getFreq(s)
% This gets the frequency for the three periodic lens modes.
clearBytes(s);
% mystr=double('PrFAaaaa');
% msg=appendCRC16(mystr);
fwrite(s,[80   114    70    65    97    97    97    97   195   200]);
reply=fread(s,11);
freq=bin2dec ([dec2bin(reply(6),8) dec2bin(reply(7),8)])/1000;
end

function setHigh(s,high)
% This sets the upper bound for current swings in the three periodic lens modes.
if high>4095
    high=4095;
end
if high<0
    high=0;
end

mystr=double('PwUA');

bytes=dec2bin(high,16);
byte1=bin2dec(bytes(1:8));
byte2=bin2dec(bytes(9:16));
byte3=0;
byte4=0;

msg=appendCRC16([mystr byte1 byte2 byte3 byte4]);
fwrite(s,msg);

clearBytes(s);
end

function setLow(s,low)
% This sets the lower bound for current swings in the three periodic lens modes.

if low>4095
    low=4095;
end
if low<0
    low=0;
end

mystr=double('PwLA');

bytes=dec2bin(low,16);
byte1=bin2dec(bytes(1:8));
byte2=bin2dec(bytes(9:16));
byte3=0;
byte4=0;

msg=appendCRC16([mystr byte1 byte2 byte3 byte4]);
fwrite(s,msg);

clearBytes(s);
end

function setFreq(s,freq)
% See getFreq. Change set -> get. ALSO: Frequency has a lower bound of 0.06
% Hz, an upper bound of 50 Hz. Setting any lower frequency will default to 1
% Hz, and the upper bound was imposed to protect the lens.

if freq<0.06
    freq=0.06;
end
if freq>50
    freq=50;
end

encfreq=1000*freq;
mystr=double('PwFA');

bytes=dec2bin(encfreq,32);
byte1=bin2dec(bytes(1:8));
byte2=bin2dec(bytes(9:16));
byte3=bin2dec(bytes(17:24));
byte4=bin2dec(bytes(25:32));

msg=appendCRC16([mystr byte1 byte2 byte3 byte4]);
fwrite(s,msg);

clearBytes(s);
end

function [max, lower, upper]=getCurrentBounds(s)
% This requests the maximum current (imposed by the hardware-- DO NOT CHANGE
% THIS), and the lower and upper current bounds (imposed by the software--
% currently maximized; really no reason to make it any smaller). 

% maxstr=double('CrMAaa');
% maxdecinput=double(maxstr);
% maxmsg=appendCRC16(maxdecinput);
% lowstr='CrLAaa';
% lowdecinput=double(lowstr);
% lowmsg=appendCRC16(lowdecinput);
% uppstr='CrUAaa';
% uppdecinput=double(uppstr);
% uppmsg=appendCRC16(uppdecinput);

clearBytes(s)
fwrite(s,[67   114    77    65    97    97   153   248]);
maxreply=fread(s,9);
max=bin2dec ([dec2bin(maxreply(4),8) dec2bin(maxreply(5),8)])/100;

clearBytes(s)
fwrite(s,[67   114    76    65    97    97   152     4]);
lowreply=fread(s,9);
lower=max*bin2dec ([dec2bin(lowreply(4),8) dec2bin(lowreply(5),8)])/4095;

clearBytes(s)
fwrite(s,[67   114    85    65    97    97   159    88]);
uppreply=fread(s,9);
upper=max*bin2dec ([dec2bin(uppreply(4),8) dec2bin(uppreply(5),8)])/4095;
end

% CURRENT READ/WRITE COMMANDS: Change the current.
function curr=getCurr(s)
% Gets the current (reported as a 12 bit number 0-4096).

clearBytes(s)
fwrite(s,[65   114     0     0   180    39]);
reply=fread(s,7);
curr=bin2dec([dec2bin(reply(2),8) dec2bin(reply(3),8)]);
end

function setCurr(s,curr)
% Sets the current. Curr is a 12 bit number, NOT in mA.

mystr=double('Aw');

if curr<0
    curr=0;
end
if curr>4095
    curr=4095;
end

level=dec2bin(curr,16);
level1=bin2dec(level(1:8));
level2=bin2dec(level(9:16));
appended=appendCRC16([mystr level1 level2]);

fwrite(s,appended);
clearBytes(s);
end

function rampCurr(s,target)
% Ramps the current in ONE of the lenses between the present current and the
% target current. Not very smooth, but it's better than simply setting the
% current to the target. This process is done via SERIAL and is thus not
% very efficient. Use the analog control method for smoother, albeit more
% roundabout, results.
curr=getCurr(s);

currentVector=round(linspace(curr,target,10));

for ii=1:length(currentVector)
    setCurr(s,currentVector(ii));
end
end

function faroCurr(s_A,s_B,atarget,btarget,res)
% Same thing as rampCurr but this is designed to ramp the currents in both
% lenses simultaneously (Faro shuffle: Google it).

clearBytes(s_A);
clearBytes(s_B);

acurr=getCurr(s_A);
bcurr=getCurr(s_B);
aVector=round(linspace(acurr,atarget,res));
bVector=round(linspace(bcurr,btarget,res));

for ii=1:res
    setCurr(s_A,aVector(ii));
    setCurr(s_B,bVector(ii));
end

end

function resetCurr(s_A,s_B)
% When we're done and we need to semi-smoothly stop the current in the lenses.

modeDC(s_A)
modeDC(s_B)

acurr=getCurr(s_A);
bcurr=getCurr(s_B);

if acurr ~= 0
    rampCurr(s_A,0)
end
if bcurr ~= 0
    rampCurr(s_B,0)
end
end

% MISCELLANEOUS / INDEV COMMANDS
function spotSizeFit(cal,imin,imax)

[M,~]=min(cal);
[~,uplimit]=min(M);
dim=size(cal);
currVector=linspace(imin,imax,dim(1));

% while 1
%     if I(uplimit)<I(uplimit+1) && I(uplimit)<3
%         break
%     end
%     uplimit=uplimit+1;
% end

M=transpose(M(1:uplimit));
xcal=transpose(currVector(1:uplimit));

fitobject=fit(xcal,M,'poly1');

figure
scatter(xcal,M)
hold on
grid on
plot(fitobject);
legend off
title('Mean Focused Spot Size Fit')
xlabel('I_B (d-mA)')
ylabel('Beam Waist (mm)')
end
function mapFit(s_Test,s_Check,AO)
modeAnalog(s_Test)
modeDC(s_Check)
clearBytes(s_Test)
clearBytes(s_Check)

analogVector=linspace(0,5,100);
currVector=linspace(0,0,100);
analogCheckVector=linspace(0,0,100);

for ii=1:100
    putsample(AO,[analogVector(ii) analogVector(ii)]);
    pause(0.1)
    currVector(ii)=getCurr(s_Test);
    analogCheckVector(ii)=getAnalog(s_Check);
    disp(ii)
end

putsample(AO,[0 0])

figure
scatter(analogVector,292.77*currVector/4095)
grid on
xlabel('Analog Input (V)')
ylabel('Lens Current (mA)')
title('Current Response in Analog Control Mode')

keyboard
end
function testfunct(s_A,s_B)
[vidobj,~]=initializeGuppyVimba(2,100)

modeDC(s_A)
modeDC(s_B)

rampCurr(s_A,140);
rampCurr(s_B,2860);
pause(0.05)
trigger(vidobj);
data1=getdata(vidobj);
pause(0.05)

rampCurr(s_A,2860);
rampCurr(s_B,140);
pause(0.05)
trigger(vidobj);
data2=getdata(vidobj);

resetCurr(s_A,s_B)

keyboard 
end
%   I. TEMPERATURE COMPENSATION AND VERIFICATION: PREPARING FOR TEMPERATURE
%   COMPENSATION

%     1. Steady state temperature data collection.
function tempdata=ssTemp(AO,s_A,s_B,imin,imax,csteps,tsteps)
% This is not compensation. This is for setting temp.
clearBytes(s_A);
clearBytes(s_B);
basename = input('ssTemp Name: ','s');
tempdata=zeros(csteps,tsteps,2);
currVector=round(linspace(imin,imax,csteps));
time=linspace(0,5*(tsteps-1),tsteps);

rampAnalog(AO,imin,imin)
pause(60)

indx=1;
for ii=1:csteps
    rampAnalog(AO,currVector(ii),currVector(ii));
    for jj=1:tsteps
        tempdata(ii,jj,1)=getTemp(s_A);
        tempdata(ii,jj,2)=getTemp(s_B);
        disp([tempdata(ii,jj,1) tempdata(ii,jj,2)])
        disp(['Progress: ' num2str (100*indx/(tsteps*csteps)) '%'])
        indx=indx+1;
        pause(5)
    end
end
rampAnalog(AO,0,0)

str=[date '\' basename '_' num2str (imin) '_' num2str (imax) '_' num2str (csteps) '_' num2str(tsteps)];
mypath = ['C:\Users\Bolton\Documents\MATLAB\Optotune\r uns\' str];
mkdir(mypath);
csvwrite(fullfile(mypath,'atempdata.csv'),tempdata(:,:,1));
csvwrite(fullfile(mypath,'btempdata.csv'),tempdata(:,:,2));

figure
h1=scatter(time,tempdata(1,:,1)); % #ok<NASGU>
hold on
grid on
for ii=2:csteps
    scatter(time,tempdata(ii,:,1))
end
xlabel('Time (s)')
ylabel('T_A (\[Degree]C)')
title('Temerature Over Time at Various Currents')

savefig(fullfile(mypath,'sstemp_a.fig'));

figure
h2=scatter(time,tempdata(1,:,2)); % #ok<NASGU>
hold on
grid on
for ii=2:csteps
    scatter(time,tempdata(ii,:,2))
end
xlabel('Time (s)')
ylabel('T_B (\[Degree]C)')
title('Temerature Over Time at Various Currents')

savefig(fullfile(mypath,'sstemp_b.fig'));
keyboard
end

%     2. Collect drift data for each lens individually (temp/curr space)
function driftMatrix=driftCalibration(AO,s_cal,s_idle,tempdata,tmin,tmax,imin,imax,tpl,cpl,tol)
% run ssTemp to generate tempdata or load an old ssTemp data file prior to
% running this function

%% DATA COLLECTION
[vidobj,~]=initializeGuppyVimba(cpl*tpl,100);

clearBytes(s_cal);
clearBytes(s_idle);
resetCurr(s_cal,s_idle);

basename = input('Run Name: ','s');

waistMatrix=zeros(tpl,cpl,3);
tempMatrix=zeros(tpl,cpl,2);
currVector=round(linspace(imin,imax,cpl));
tempVector=round(linspace(tmin,tmax,tpl),4);
click=0;

destempfit=getDesiredTempFit(tempdata,imin,imax);

modeAnalog(s_cal)
modeDC(s_idle)

for ii=1:tpl
    for jj=1:cpl
        disp(['Progress: ' num2str (100*click/(cpl*tpl)) '%']);
        tempCal=getTemp(s_cal);
        while (tempCal < tempVector(ii) - tol) || (tempCal > tempVector(ii) + tol)
            indx=0;
            while tempCal < (tempVector(ii) - tol)
                reqcurr=desTemp(destempfit,tempVector(ii)+(1+indx/50),imin,imax);
                rampAnalog(AO,reqcurr,reqcurr)
                pause(1);
                tempCal=getTemp(s_cal);
                disp(['Temperature ' get (s_cal,'Name') ': ' num2str (tempCal) ' \[Degree]C']);
                indx=indx+1;
            end
            indx=0;
            while tempCal > (tempVector(ii) + tol)
                reqcurr=desTemp(destempfit,tempVector(ii)-(1+indx/50),imin,imax);
                rampAnalog(AO,reqcurr,reqcurr)
                pause(1);
                tempCal=getTemp(s_cal);
                disp(['Temperature ' get (s_cal,'Name') ': ' num2str (tempCal) ' \[Degree]C']);
                indx=indx+1;
            end
        end
        rampAnalog(AO,currVector(jj),currVector(jj))
        pause(0.05)
        tempCal=getTemp(s_cal);
        tempIdl=getTemp(s_idle);
        trigger(vidobj);
        maintaincurr=desTemp(destempfit,tempVector(ii),imin,imax);
        rampAnalog(AO,maintaincurr,maintaincurr)
        data=getdata(vidobj);
        [waist_x,waist_y]=guppyAnalysis_rev6(data);
        waistMatrix(ii,jj,1)=waist_x;
        waistMatrix(ii,jj,2)=waist_y;
        waistMatrix(ii,jj,3)=mean(waistMatrix(ii,jj,1:2));
        tempMatrix(ii,jj,1)=tempCal;
        tempMatrix(ii,jj,2)=tempIdl;
        
        click=click+1;
    end
end

disp('Progress: 100%');

resetCurr(s_cal,s_idle);

driftMatrix=waistMatrix(:,:,1);

%% SAVE DATA TO FILE
str=[date '\' basename '_' num2str (tmin) '_' num2str (tmax) '_' num2str (imin) '_' num2str (imax) '_' num2str (tpl) '_' num2str (cpl) '_' num2str(tol)];
mypath = ['C:\Users\Bolton\Documents\MATLAB\Optotune\r uns\' str];
mkdir(mypath);
waistname = '_waists.csv';
tempname = '_temps.csv';

csvwrite(fullfile(mypath,['x' waistname]),waistMatrix(:,:,1));
csvwrite(fullfile(mypath,['A' tempname]),tempMatrix(:,:,1));
csvwrite(fullfile(mypath,['y' waistname]),waistMatrix(:,:,2));
csvwrite(fullfile(mypath,['B' tempname]),tempMatrix(:,:,2));
csvwrite(fullfile(mypath,['avg' waistname]),waistMatrix(:,:,3));

%% Plots
figure
pcolor(transpose(waistMatrix(:,:,1)));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'XTick', linspace(1,tpl,9));
set(gca, 'YTick', linspace(1,cpl,9));
set(gca, 'XTickLabel', linspace(tmin,tmax,9));
set(gca, 'YTickLabel', round(296*linspace (imin,imax,9)/4096));
ylabel(['I_' get (s_cal,'Name') ' (mA)'])
xlabel(['T_' get (s_cal,'Name') ' (\[Degree]C)'])
caxis([min(min(waistMatrix(:,:,1))) max(max(waistMatrix(:,:,1)))])
title('Guppy X Spot Diameter (mm)')

savefig(fullfile(mypath,'xspot.fig'));

figure
pcolor(transpose(waistMatrix(:,:,2)));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'XTick', linspace(1,tpl,9));
set(gca, 'YTick', linspace(1,cpl,9));
set(gca, 'XTickLabel', linspace(tmin,tmax,9));
set(gca, 'YTickLabel', round(296*linspace (imin,imax,9)/4096));
ylabel(['I_' get (s_cal,'Name') ' (mA)'])
xlabel(['T_' get (s_cal,'Name') ' (\[Degree]C)'])
caxis([min(min(waistMatrix(:,:,2))) max(max(waistMatrix(:,:,2)))])
title('Guppy Y Spot Diameter (mm)')

savefig(fullfile(mypath,'yspot.fig'));

figure
pcolor(transpose(waistMatrix(:,:,3)));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'XTick', linspace(1,tpl,9));
set(gca, 'YTick', linspace(1,cpl,9));
set(gca, 'XTickLabel', linspace(tmin,tmax,9));
set(gca, 'YTickLabel', round(296*linspace (imin,imax,9)/4096));
ylabel(['I_' get (s_cal,'Name') ' (mA)'])
xlabel(['T_' get (s_cal,'Name') ' (\[Degree]C)'])
caxis([min(min(waistMatrix(:,:,3))) max(max(waistMatrix(:,:,3)))])
title('Guppy Avg Spot Diameter (mm)')

savefig(fullfile(mypath,'avgspot.fig'));

figure
subplot(2,2,1)
pcolor(transpose(waistMatrix(:,:,1)));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'XTick', linspace(1,tpl,6));
set(gca, 'YTick', linspace(1,cpl,4));
set(gca, 'XTickLabel', linspace(tmin,tmax,6));
set(gca, 'YTickLabel', round(296*linspace (imin,imax,4)/4096));
ylabel(['I_' get (s_cal,'Name') ' (mA)'])
xlabel(['T_' get (s_cal,'Name') ' (\[Degree]C)'])
caxis([min(min(waistMatrix(:,:,1))) max(max(waistMatrix(:,:,1)))])
title('Guppy X Spot Diameter (mm)')

subplot(2,2,2)
pcolor(transpose(waistMatrix(:,:,2)));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'XTick', linspace(1,tpl,6));
set(gca, 'YTick', linspace(1,cpl,4));
set(gca, 'XTickLabel', linspace(tmin,tmax,6));
set(gca, 'YTickLabel', round(296*linspace (imin,imax,4)/4096));
ylabel(['I_' get (s_cal,'Name') ' (mA)'])
xlabel(['T_' get (s_cal,'Name') ' (\[Degree]C)'])
caxis([min(min(waistMatrix(:,:,2))) max(max(waistMatrix(:,:,2)))])
title('Guppy Y Spot Diameter (mm)')

subplot(2,2,3)
imagesc(transpose(tempMatrix(:,:,1)));
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'XTick', linspace(1,tpl,6));
set(gca, 'YTick', linspace(1,cpl,4));
set(gca, 'XTickLabel', linspace(tmin,tmax,6));
set(gca, 'YTickLabel', round(296*linspace (imin,imax,4)/4096));
ylabel(['I_' get (s_cal,'Name') ' (mA)'])
xlabel(['T_' get (s_cal,'Name') ' (\[Degree]C)'])
caxis([min(min(tempMatrix(:,:,1))) max(max(tempMatrix(:,:,1)))])
title(['Temperature of Lens ' get (s_cal,'Name') ' (\[Degree]C)'])

subplot(2,2,4)
imagesc(transpose(tempMatrix(:,:,2)));
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'XTick', linspace(1,tpl,6));
set(gca, 'YTick', linspace(1,cpl,4));
set(gca, 'XTickLabel', linspace(tmin,tmax,6));
set(gca, 'YTickLabel', round(296*linspace (imin,imax,4)/4096));
ylabel(['I_' get (s_cal,'Name') ' (mA)'])
xlabel(['T_' get (s_cal,'Name') ' (\[Degree]C)'])
caxis([min(min(tempMatrix(:,:,2))) max(max(tempMatrix(:,:,2)))])
title(['Temperature of Lens ' get (s_idle,'Name') ' (\[Degree]C)'])

savefig(fullfile(mypath,'2by2.fig'));

keyboard
end
function destempfit=getDesiredTempFit(tempdata,imin,imax)
% This is not compensation. This is for setting temp.
[csteps,~]=size(tempdata);
avgt=zeros(csteps,1);
currVector=transpose(round(linspace(imin,imax,csteps)));

for ii=1:csteps
    avgt(ii)=max(tempdata(ii,:));
end

destempfit=fit(avgt,currVector,'smoothingspline');
end
function reqcurr=desTemp(destempfit,destemp,imin,imax)
% This is not compensation. This is for setting temp.
reqcurr=round(destempfit(destemp));
if reqcurr < imin
    reqcurr=imin;
end
if reqcurr > imax
    reqcurr=imax;
end
end

%     3. Verify that the current interpolation method (compCurrInterp) can
%       remove temperature dependence from optical behavior.
function driftTestMatrix=driftCalTest(AO,s_cal,s_idle,tempdata,driftMatrix,tmin,tmax,imin,imax,tpl,cpl,tol)
%% DATA COLLECTION
[vidobj,~]=initializeGuppyVimba(cpl*tpl,100);

clearBytes(s_cal);
clearBytes(s_idle);
resetCurr(s_cal,s_idle);

basename = input('Run Name: ','s');

waistMatrix=zeros(tpl,cpl,3);
tempMatrix=zeros(tpl,cpl,2);
currVector=round(linspace(imin,imax,cpl));
tempVector=round(linspace(tmin,tmax,tpl),4);
click=0;

destempfit=getDesiredTempFit(tempdata,imin,imax);

modeAnalog(s_cal)
modeDC(s_idle)

for ii=1:tpl
    for jj=1:cpl
        disp(['Progress: ' num2str (100*click/(cpl*tpl)) '%']);
        tempCal=getTemp(s_cal);
        while (tempCal < tempVector(ii) - tol) || (tempCal > tempVector(ii) + tol)
            indx=0;
            while tempCal < tempVector(ii) - tol
                reqcurr=desTemp(destempfit,tempVector(ii)+(1+indx/50),imin,imax);
                rampAnalog(AO,reqcurr,reqcurr);
                pause(1);
                tempCal=getTemp(s_cal);
                disp(['Temperature ' get (s_cal,'Name') ': ' num2str (tempCal) ' \[Degree]C']);
                indx=indx+1;
            end
            indx=0;
            while tempCal > tempVector(ii) + tol
                reqcurr=desTemp(destempfit,tempVector(ii)-(1+indx/50),imin,imax);
                rampAnalog(AO,reqcurr,reqcurr);
                pause(1);
                tempCal=getTemp(s_cal);
                disp(['Temperature ' get (s_cal,'Name') ': ' num2str (tempCal) ' \[Degree]C']);
                indx=indx+1;
            end
        end
        outputcurr=compCurrInterp(s_cal,driftMatrix,tmin,tmax,imin,imax,currVector(jj));
        rampAnalog(AO,outputcurr,outputcurr)
        pause(0.05)
        tempCal=getTemp(s_cal);
        tempIdl=getTemp(s_idle);
        trigger(vidobj);
        maintaincurr=desTemp(destempfit,tempVector(ii),imin,imax);
        rampAnalog(AO,maintaincurr,maintaincurr);
        data=getdata(vidobj);
        [waist_x,waist_y]=guppyAnalysis_rev6(data);
        waistMatrix(ii,jj,1)=waist_x;
        waistMatrix(ii,jj,2)=waist_y;
        waistMatrix(ii,jj,3)=mean(waistMatrix(ii,jj,1:2));
        tempMatrix(ii,jj,1)=tempCal;
        tempMatrix(ii,jj,2)=tempIdl;
        click=click+1;
    end
end

disp('Progress: 100%');

resetCurr(s_cal,s_idle);

driftTestMatrix=waistMatrix(:,:,1);

%% SAVE DATA TO FILE
str=[date '\' basename '_' num2str (tmin) '_' num2str (tmax) '_' num2str (imin) '_' num2str (imax) '_' num2str (tpl) '_' num2str (cpl) '_' num2str(tol)];
mypath = ['C:\Users\Bolton\Documents\MATLAB\Optotune\r uns\' str];
mkdir(mypath);
waistname = '_waists.csv';
tempname = '_temps.csv';

csvwrite(fullfile(mypath,['x' waistname]),waistMatrix(:,:,1));
csvwrite(fullfile(mypath,['A' tempname]),tempMatrix(:,:,1));
csvwrite(fullfile(mypath,['y' waistname]),waistMatrix(:,:,2));
csvwrite(fullfile(mypath,['B' tempname]),tempMatrix(:,:,2));
csvwrite(fullfile(mypath,['avg' waistname]),waistMatrix(:,:,3));

%% Plots
figure
pcolor(transpose(waistMatrix(:,:,1)));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'XTick', linspace(1,tpl,9));
set(gca, 'YTick', linspace(1,cpl,9));
set(gca, 'XTickLabel', linspace(tmin,tmax,9));
set(gca, 'YTickLabel', round(296*linspace (imin,imax,9)/4096));
ylabel('I_A (mA)')
xlabel('T_A (\[Degree]C)')
caxis([min(min(waistMatrix(:,:,1))) max(max(waistMatrix(:,:,1)))])
title('Guppy X Spot Diameter (mm)')

savefig(fullfile(mypath,'xspot.fig'));

figure
pcolor(transpose(waistMatrix(:,:,2)));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'XTick', linspace(1,tpl,9));
set(gca, 'YTick', linspace(1,cpl,9));
set(gca, 'XTickLabel', linspace(tmin,tmax,9));
set(gca, 'YTickLabel', round(296*linspace (imin,imax,9)/4096));
ylabel('I_A (mA)')
xlabel('T_A (\[Degree]C)')
caxis([min(min(waistMatrix(:,:,2))) max(max(waistMatrix(:,:,2)))])
title('Guppy Y Spot Diameter (mm)')

savefig(fullfile(mypath,'yspot.fig'));

figure
pcolor(transpose(waistMatrix(:,:,3)));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'XTick', linspace(1,tpl,9));
set(gca, 'YTick', linspace(1,cpl,9));
set(gca, 'XTickLabel', linspace(tmin,tmax,9));
set(gca, 'YTickLabel', round(296*linspace (imin,imax,9)/4096));
ylabel('I_A (mA)')
xlabel('T_A (\[Degree]C)')
caxis([min(min(waistMatrix(:,:,3))) max(max(waistMatrix(:,:,3)))])
title('Guppy Avg Spot Diameter (mm)')

savefig(fullfile(mypath,'avgspot.fig'));

figure
subplot(2,2,1)
pcolor(transpose(waistMatrix(:,:,1)));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'XTick', linspace(1,tpl,6));
set(gca, 'YTick', linspace(1,cpl,4));
set(gca, 'XTickLabel', linspace(tmin,tmax,6));
set(gca, 'YTickLabel', round(296*linspace (imin,imax,4)/4096));
ylabel('I_A (mA)')
xlabel('T_A (\[Degree]C)')
caxis([min(min(waistMatrix(:,:,1))) max(max(waistMatrix(:,:,1)))])
title('Guppy X Spot Diameter (mm)')

subplot(2,2,2)
pcolor(transpose(waistMatrix(:,:,2)));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'XTick', linspace(1,tpl,6));
set(gca, 'YTick', linspace(1,cpl,4));
set(gca, 'XTickLabel', linspace(tmin,tmax,6));
set(gca, 'YTickLabel', round(296*linspace (imin,imax,4)/4096));
ylabel('I_A (mA)')
xlabel('T_A (\[Degree]C)')
caxis([min(min(waistMatrix(:,:,2))) max(max(waistMatrix(:,:,2)))])
title('Guppy Y Spot Diameter (mm)')

subplot(2,2,3)
imagesc(transpose(tempMatrix(:,:,1)));
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'XTick', linspace(1,tpl,6));
set(gca, 'YTick', linspace(1,cpl,4));
set(gca, 'XTickLabel', linspace(tmin,tmax,6));
set(gca, 'YTickLabel', round(296*linspace (imin,imax,4)/4096));
ylabel('I_A (mA)')
xlabel('T_A (\[Degree]C)')
caxis([min(min(tempMatrix(:,:,1))) max(max(tempMatrix(:,:,1)))])
title('Temperature of Lens A (\[Degree]C)')

subplot(2,2,4)
imagesc(transpose(tempMatrix(:,:,2)));
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'XTick', linspace(1,tpl,6));
set(gca, 'YTick', linspace(1,cpl,4));
set(gca, 'XTickLabel', linspace(tmin,tmax,6));
set(gca, 'YTickLabel', round(296*linspace (imin,imax,4)/4096));
ylabel('I_A (mA)')
xlabel('T_A (\[Degree]C)')
caxis([min(min(tempMatrix(:,:,2))) max(max(tempMatrix(:,:,2)))])
title('Temperature of Lens B (\[Degree]C)')

savefig(fullfile(mypath,'2by2.fig'));

stop(vidobj)
keyboard
end
function outputcurr=compCurrInterp(s,driftMatrix,tmin,tmax,imin,imax,curr)
clearBytes(s)
temp=getTemp(s);

if temp<tmin
    temp=tmin;
end
if temp>tmax
    temp=tmax;
end

if curr<imin
    curr=imin;
end
if curr>imax
    curr=imax;
end

basetemp=26.5;

dataMatrix=transpose(driftMatrix);
[cpl, tpl]=size(dataMatrix);

tempVector=transpose(linspace(tmin,tmax,tpl));
currVector=transpose(linspace(imin,imax,cpl));
[~,I]=min(dataMatrix);
mindexfit=fit(tempVector,transpose(I),'smoothingspline');

minindx=round(mindexfit(temp));
deltac=curr-currVector(minindx);

[tempMatrix,currMatrix]=meshgrid(tempVector,currVector);
desiredspot=interp2(tempMatrix,currMatrix,dataMatrix,basetemp,curr);

if deltac<=0
    isoCurrVector=currVector(1:minindx);
    isoSpotVector=zeros(minindx,1);
    for ii=1:minindx
        isoSpotVector(ii)=interp2(tempMatrix,currMatrix,dataMatrix,temp,isoCurrVector(ii));
    end
    if desiredspot < min(isoSpotVector)
        desiredspot = min(isoSpotVector);
    end
    if desiredspot > max(isoSpotVector)
        desiredspot = max(isoSpotVector);
    end
    outputcurr=interp1(isoSpotVector,isoCurrVector,desiredspot);
end

if deltac>0
    isoCurrVector=currVector(minindx:cpl);
    isoSpotVector=zeros(cpl-minindx+1,1);
    for ii=1:cpl-minindx+1
        isoSpotVector(ii)=interp2(tempMatrix,currMatrix,dataMatrix,temp,isoCurrVector(ii));
    end
    if desiredspot < min(isoSpotVector)
        desiredspot = min(isoSpotVector);
    end
    if desiredspot > max(isoSpotVector)
        desiredspot = max(isoSpotVector);
    end
    outputcurr=interp1(isoSpotVector,isoCurrVector,desiredspot);
end
end
function driftTest=driftCalHold(AO,s_cal,s_idle,tempdata,driftMatrix,tmin,tmax,imin,imax)
end

%   II. ACQUIRE (UN)COMPENSATED CALIBRATION DATA IN CURRENT-CURRENT SPACE
%     
%     1. Randomly collect spot size data in the I_A-I_B current space.
function cal=getCalAnalog(AO,s_A,s_B,~,tmax,imin,imax,ppl)
%% SETUP
clearBytes(s_A);
clearBytes(s_B);

atarget=0;
btarget=0;

[vidobj,~]=initializeGuppyVimba(ppl*ppl,100);

basename = input('Run Name: ','s');

modeAnalog(s_A);
modeAnalog(s_B);

rampAnalog(AO,imin,imin);

currVector=round(linspace(imin, imax, ppl)); % imax is given on a a scale from 0 to 4095
waistMatrix=zeros(ppl,ppl,3);
tempMatrix=zeros(ppl,ppl,2);

randindx=transpose(randperm(ppl*ppl));
refMatrix=zeros(ppl*ppl,2);

indx=1;
for ii=1:ppl
    for jj=1:ppl
        refMatrix(indx,1)=ii;
        refMatrix(indx,2)=jj;
        indx=indx+1;
    end
end

%% WARMUP
tempA=getTemp(s_A);
tempB=getTemp(s_B);

while (tempA < 23.8460) || (tempB < 22.8858)
    if tempA < 23.8460
        setA=4000;
    else
        setA=2000;
    end
    if tempB < 22.8858
        setB=4000;
    else
        setB=2000;
    end
    rampAnalog(AO,setA,setB)
    pause(1)
    tempA=getTemp(s_A);
    tempB=getTemp(s_B);
    disp(['Temperature A: ' num2str (tempA) ' \[Degree]C']);
    disp(['Temperature B: ' num2str (tempB) ' \[Degree]C']);
end

rampAnalog(AO,imin,imin)

%% DATA COLLECTION
tic
indx=1;
for ii=1:ppl
    for jj=1:ppl
        tempA=getTemp(s_A);
        tempB=getTemp(s_B);
        disp(['Progress: ' num2str (100*(indx)/(ppl)^2) '%']);
        time=toc;
        disp(['Time: ' num2str (time) ' s']);
        eta=((ppl)^2-(ppl*(ii-1)+jj))*time/(ppl*(ii-1)+jj);
        [q, sec] = quorem(round(sym(eta)), sym(60));
        [hr, minutes]=quorem(sym(q), sym(60));
        disp(['Est:  ' num2str (int8(hr)) ':' num2str (int8(minutes)) ':' num2str(int8(sec))]);
        disp(' ');
        if (tempA >= tmax) || (tempB >=tmax)
            rampAnalog(imin,imin)
            while (tempA >= tmax) || (tempB >=tmax)
                pause(5)
                tempA=getTemp(s_A);
                tempB=getTemp(s_B);
                disp(['Temperature A: ' num2str (tempA) ' \[Degree]C']);
                disp(['Temperature B: ' num2str (tempB) ' \[Degree]C']);
            end
        end
        indxA=refMatrix(randindx(ppl*(ii-1)+jj),1);
        indxB=refMatrix(randindx(ppl*(ii-1)+jj),2);
        rampAnalog(AO,currVector(indxA),currVector(indxB));
        pause(0.05)
        tempA=getTemp(s_A);
        tempB=getTemp(s_B);
        trigger(vidobj);
        rampAnalog(AO,imin,imin);
        data=getdata(vidobj);
        [waist_x,waist_y]=guppyAnalysis_rev6(data);
%         if waist_x > 1 || waist_y > 1
%             keyboard
%         end
        waistMatrix(indxA,indxB,1)=waist_x;
        waistMatrix(indxA,indxB,2)=waist_y;
        waistMatrix(indxA,indxB,3)=mean([waist_x waist_y]);
        tempMatrix(indxA,indxB,1)=tempA;
        tempMatrix(indxA,indxB,2)=tempB;
        indx=indx+1;
    end
end

rampAnalog(AO,0,0)
resetCurr(s_A,s_B);
cal=waistMatrix(:,:,3);

[M,I]=min(waistMatrix(:,:,2));

uplimit=1;

while 1
    if I(uplimit)<I(uplimit+1) && I(uplimit)<3
        break
    end
    uplimit=uplimit+1;
end

M=transpose(M(1:uplimit));
xcal=transpose(currVector(1:uplimit));

%% SAVE DATA TO FILES
str=[date '\' basename '_' num2str (imin) '_' num2str (imax) '_' num2str(ppl)];
mypath = ['C:\Users\Bolton\Documents\MATLAB\Optotune\r uns\' str];
mkdir(mypath);
waistname = '_waists.csv';
tempname = '_temps.csv';

csvwrite(fullfile(mypath,['avg' waistname]),waistMatrix(:,:,3));
csvwrite(fullfile(mypath,['x' waistname]),waistMatrix(:,:,1));
csvwrite(fullfile(mypath,['A' tempname]),tempMatrix(:,:,1));
csvwrite(fullfile(mypath,['y' waistname]),waistMatrix(:,:,2));
csvwrite(fullfile(mypath,['B' tempname]),tempMatrix(:,:,2));

%% PLOTS
figure
pcolor(waistMatrix(:,:,1));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'YTick', linspace(1,ppl,10));
set(gca, 'XTick', linspace(1,ppl,10));
set(gca, 'YTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
set(gca, 'XTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
xlabel('I_B (mA)')
ylabel('I_A (mA)')
title('Guppy X Spot Diameter (mm)')

savefig(fullfile(mypath,'xspot.fig'));

figure
pcolor(waistMatrix(:,:,2));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'YTick', linspace(1,ppl,10));
set(gca, 'XTick', linspace(1,ppl,10));
set(gca, 'YTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
set(gca, 'XTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
xlabel('I_B (mA)')
ylabel('I_A (mA)')
title('Guppy Y Spot Diameter (mm)')

savefig(fullfile(mypath,'yspot.fig'));

figure
pcolor(waistMatrix(:,:,3));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'YTick', linspace(1,ppl,10));
set(gca, 'XTick', linspace(1,ppl,10));
set(gca, 'YTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
set(gca, 'XTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
xlabel('I_B (mA)')
ylabel('I_A (mA)')
% caxis([0 1.5])
title('Mean of X and Y Guppy Spot Diameters (mm)')

savefig(fullfile(mypath,'avgspot.fig'));

figure
subplot(2,2,1)
pcolor(waistMatrix(:,:,1));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'YTick', linspace(1,ppl,10));
set(gca, 'XTick', linspace(1,ppl,10));
set(gca, 'YTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
set(gca, 'XTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
xlabel('I_B (mA)')
ylabel('I_A (mA)')
% caxis([0 1.5])
title('Guppy X Spot Diameter (mm)')

subplot(2,2,2)
pcolor(waistMatrix(:,:,2));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'YTick', linspace(1,ppl,10));
set(gca, 'XTick', linspace(1,ppl,10));
set(gca, 'YTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
set(gca, 'XTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
xlabel('I_B (mA)')
ylabel('I_A (mA)')
% caxis([0 1.5])
title('Guppy Y Spot Diameter (mm)')

subplot(2,2,3)
imagesc(tempMatrix(:,:,1));
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'YTick', linspace(1,ppl,10));
set(gca, 'XTick', linspace(1,ppl,10));
set(gca, 'YTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
set(gca, 'XTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
xlabel('I_B (mA)')
ylabel('I_A (mA)')
caxis([min(min(tempMatrix(:,:,1))) max(max(tempMatrix(:,:,1)))])
title('Temperature of Lens A (\[Degree]C)')

subplot(2,2,4)
imagesc(tempMatrix(:,:,2));
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottomdele
set(gca, 'YTick', linspace(1,ppl,10));
set(gca, 'XTick', linspace(1,ppl,10));
set(gca, 'YTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
set(gca, 'XTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
xlabel('I_B (mA)')
ylabel('I_A (mA)')
caxis([min(min(tempMatrix(:,:,2))) max(max(tempMatrix(:,:,2)))])
title('Temperature of Lens B (\[Degree]C)')

savefig(fullfile(mypath,'2by2.fig'));

figure
scatter(xcal,M);
hold on
grid on
xlabel('I_B (d-mA)');
ylabel('Beam Waist (mm)');
title('Uncompensated Calibration: Average Focused Spot Size')

savefig(fullfile(mypath,'avgfocus.fig'));

keyboard
end
function compcal=getCompCalAnalog(AO,s_A,s_B,tmin,tmax,imin,imax,ppl,adrift,bdrift)
%% SETUP
clearBytes(s_A);
clearBytes(s_B);

[vidobj,~]=initializeGuppyVimba(ppl*ppl,100);

basename = input('Run Name: ','s');

modeAnalog(s_A);
modeAnalog(s_B);

rampAnalog(AO,imin,imin)

currVector=round(linspace(imin, imax, ppl)); % imax is given on a a scale from 0 to 4095
waistMatrix=zeros(ppl,ppl,3);
tempMatrix=zeros(ppl,ppl,2);

randindx=transpose(randperm(ppl*ppl));
refMatrix=zeros(ppl*ppl,2);

indx=1;
for ii=1:ppl
    for jj=1:ppl
        refMatrix(indx,1)=ii;
        refMatrix(indx,2)=jj;
        indx=indx+1;
    end
end

%% WARMUP
tempA=getTemp(s_A);
tempB=getTemp(s_B);

while (tempA < 23.8460) || (tempB < 22.8858)
    if tempA < 23.8460
        setA=4000;
    else
        setA=2000;
    end
    if tempB < 22.8858
        setB=4000;
    else
        setB=2000;
    end
    rampAnalog(AO,setA,setB)
    pause(1)
    tempA=getTemp(s_A);
    tempB=getTemp(s_B);
    disp(['Temperature A: ' num2str (tempA) ' \[Degree]C']);
    disp(['Temperature B: ' num2str (tempB) ' \[Degree]C']);
end
rampAnalog(AO,imin,imin)

%% DATA COLLECTION
tic
indx=1;
for ii=1:ppl
    for jj=1:ppl
        tempA=getTemp(s_A);
        tempB=getTemp(s_B);
        disp(['Progress: ' num2str (100*(indx)/(ppl)^2) '%']);
        time=toc;
        disp(['Time: ' num2str (time) ' s']);
        eta=((ppl)^2-(ppl*(ii-1)+jj))*time/(ppl*(ii-1)+jj);
        [q, sec] = quorem(round(sym(eta)), sym(60));
        [hr, minutes]=quorem(sym(q), sym(60));
        disp(['Est:  ' num2str (int8(hr)) ':' num2str (int8(minutes)) ':' num2str(int8(sec))]);
        disp(' ');
        if (tempA >= tmax) || (tempB >=tmax)
            rampAnalog(AO,140,140)
            while (tempA >= tmax) || (tempB >=tmax)
                pause(5);
                tempA=getTemp(s_A);
                tempB=getTemp(s_B);
                disp(['Temperature A: ' num2str (tempA) ' \[Degree]C']);
                disp(['Temperature B: ' num2str (tempB) ' \[Degree]C']);
            end
        end
        indxA=refMatrix(randindx(ppl*(ii-1)+jj),1);
        indxB=refMatrix(randindx(ppl*(ii-1)+jj),2);
        currA=compCurrInterp(s_A,adrift,tmin,tmax,imin,imax,currVector(indxA));
        currB=compCurrInterp(s_B,bdrift,tmin,tmax,imin,imax,currVector(indxB));
        rampAnalog(AO,currA,currVector(indxB));
        pause(0.05)
        tempA=getTemp(s_A);
        tempB=getTemp(s_B);
        trigger(vidobj);
        rampAnalog(AO,imin,imin);
        data=getdata(vidobj);
        [waist_x,waist_y]=guppyAnalysis_rev6(data);
        waistMatrix(indxA,indxB,1)=waist_x;
        waistMatrix(indxA,indxB,2)=waist_y;
        waistMatrix(indxA,indxB,3)=mean([waist_x waist_y]);
        tempMatrix(indxA,indxB,1)=tempA;
        tempMatrix(indxA,indxB,2)=tempB;
        indx=indx+1;
    end
end

resetCurr(s_A,s_B);
compcal=waistMatrix(:,:,3);

[M,I]=min(waistMatrix(:,:,3));

uplimit=1;

while 1
    if I(uplimit)<I(uplimit+1) && I(uplimit)<3
        break
    end
    uplimit=uplimit+1;
end

M=transpose(M(1:uplimit));
xcal=transpose(currVector(1:uplimit));

%% SAVE DATA TO FILES
str=[date '\' basename '_' num2str (imin) '_' num2str (imax) '_' num2str(ppl)];
mypath = ['C:\Users\Bolton\Documents\MATLAB\Optotune\r uns\' str];
mkdir(mypath);
waistname = '_waists.csv';
tempname = '_temps.csv';

csvwrite(fullfile(mypath,['avg' waistname]),waistMatrix(:,:,3));
csvwrite(fullfile(mypath,['x' waistname]),waistMatrix(:,:,1));
csvwrite(fullfile(mypath,['A' tempname]),tempMatrix(:,:,1));
csvwrite(fullfile(mypath,['y' waistname]),waistMatrix(:,:,2));
csvwrite(fullfile(mypath,['B' tempname]),tempMatrix(:,:,2));

%% PLOTS
figure
pcolor(waistMatrix(:,:,1));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'YTick', linspace(1,ppl,10));
set(gca, 'XTick', linspace(1,ppl,10));
set(gca, 'YTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
set(gca, 'XTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
xlabel('I_B (mA)')
ylabel('I_A (mA)')
% caxis([min(min(waistMatrix(:,:,1))) max(max(waistMatrix(:,:,1)))])
title('Guppy X Spot Diameter (mm)')

savefig(fullfile(mypath,'xspot.fig'));

figure
pcolor(waistMatrix(:,:,2));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'YTick', linspace(1,ppl,10));
set(gca, 'XTick', linspace(1,ppl,10));
set(gca, 'YTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
set(gca, 'XTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
xlabel('I_B (mA)')
ylabel('I_A (mA)')
% caxis([min(min(waistMatrix(:,:,1))) max(max(waistMatrix(:,:,1)))])
title('Guppy Y Spot Diameter (mm)')

savefig(fullfile(mypath,'yspot.fig'));

figure
pcolor(waistMatrix(:,:,3));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'YTick', linspace(1,ppl,10));
set(gca, 'XTick', linspace(1,ppl,10));
set(gca, 'YTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
set(gca, 'XTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
xlabel('I_B (mA)')
ylabel('I_A (mA)')
% caxis([0 1.5])
title('Mean of X and Y Guppy Spot Diameters (mm)')

savefig(fullfile(mypath,'meanspot.fig'));

figure
subplot(2,2,1)
pcolor(waistMatrix(:,:,1));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'YTick', linspace(1,ppl,10));
set(gca, 'XTick', linspace(1,ppl,10));
set(gca, 'YTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
set(gca, 'XTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
xlabel('I_B (mA)')
ylabel('I_A (mA)')
% caxis([0 1.5])
title('Guppy X Spot Diameter (mm)')

subplot(2,2,2)
pcolor(waistMatrix(:,:,2));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'YTick', linspace(1,ppl,10));
set(gca, 'XTick', linspace(1,ppl,10));
set(gca, 'YTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
set(gca, 'XTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
xlabel('I_B (mA)')
ylabel('I_A (mA)')
% caxis([0 1.5])
title('Guppy Y Spot Diameter (mm)')

subplot(2,2,3)
imagesc(tempMatrix(:,:,1));
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'YTick', linspace(1,ppl,10));
set(gca, 'XTick', linspace(1,ppl,10));
set(gca, 'YTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
set(gca, 'XTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
xlabel('I_B (mA)')
ylabel('I_A (mA)')
caxis([min(min(tempMatrix(:,:,1))) max(max(tempMatrix(:,:,1)))])
title('Temperature of Lens A (\[Degree]C)')

subplot(2,2,4)
imagesc(tempMatrix(:,:,2));
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'YTick', linspace(1,ppl,10));
set(gca, 'XTick', linspace(1,ppl,10));
set(gca, 'YTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
set(gca, 'XTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
xlabel('I_B (mA)')
ylabel('I_A (mA)')
caxis([min(min(tempMatrix(:,:,2))) max(max(tempMatrix(:,:,2)))])
title('Temperature of Lens B (\[Degree]C)')

savefig(fullfile(mypath,'2by2.fig'));

figure
scatter(xcal,M);
hold on
grid on
xlabel('I_B (d-mA)');
ylabel('Beam Waist (mm)');
title('Compensated Calibration: Average Focused Spot Size')

savefig(fullfile(mypath,'avgfocus.fig'));

keyboard
end 
function cal=getCalRand(s_A,s_B,~,tmax,imin,imax,ppl)
%% SETUP
clearBytes(s_A);
clearBytes(s_B);

[vidobj,~]=initializeGuppy(ppl*ppl);

basename = input('Run Name: ','s');

rampCurr(s_A,imin);
rampCurr(s_B,imin);

currVector=round(linspace(imin, imax, ppl)); % imax is given on a a scale from 0 to 4095
waistMatrix=zeros(ppl,ppl,3);
tempMatrix=zeros(ppl,ppl,2);

randindx=transpose(randperm(ppl*ppl));
refMatrix=zeros(ppl*ppl,2);

indx=1;
for ii=1:ppl
    for jj=1:ppl
        refMatrix(indx,1)=ii;
        refMatrix(indx,2)=jj;
        indx=indx+1;
    end
end

%% WARMUP
tempA=getTemp(s_A);
tempB=getTemp(s_B);

while (tempA < 27) || (tempB < 26.5)
    if tempA < 27
        rampCurr(s_A,4000);
    else
        rampCurr(s_A,2000);
    end
    if tempB < 26.5
        rampCurr(s_B,4000);
    else
        rampCurr(s_B,2000);
    end
    pause(1)
    tempA=getTemp(s_A);
    tempB=getTemp(s_B);
    disp(['Temperature A: ' num2str (tempA) ' \[Degree]C']);
    disp(['Temperature B: ' num2str (tempB) ' \[Degree]C']);
end
rampCurr(s_A,imin);
rampCurr(s_B,imin);

%% DATA COLLECTION
tic
indx=1;
for ii=1:ppl
    for jj=1:ppl
        tempA=getTemp(s_A);
        tempB=getTemp(s_B);
        disp(['Progress: ' num2str (100*(indx)/(ppl)^2) '%']);
        time=toc;
        disp(['Time: ' num2str (time) ' s']);
        eta=((ppl)^2-(ppl*(ii-1)+jj))*time/(ppl*(ii-1)+jj);
        [q, sec] = quorem(round(sym(eta)), sym(60));
        [hr, minutes]=quorem(sym(q), sym(60));
        disp(['Est:  ' num2str (int8(hr)) ':' num2str (int8(minutes)) ':' num2str(int8(sec))]);
        disp(' ');
        if (tempA >= tmax) || (tempB >=tmax)
            rampCurr(s_A,140);
            rampCurr(s_B,140);
            while (tempA >= tmax) || (tempB >=tmax)
                pause(5)
                tempA=getTemp(s_A);
                tempB=getTemp(s_B);
                disp(['Temperature A: ' num2str (tempA) ' \[Degree]C']);
                disp(['Temperature B: ' num2str (tempB) ' \[Degree]C']);
            end
        end
        indxA=refMatrix(randindx(ppl*(ii-1)+jj),1);
        indxB=refMatrix(randindx(ppl*(ii-1)+jj),2);
        rampCurr(s_A,currVector(indxA));
        rampCurr(s_B,currVector(indxB));
        pause(0.05)
        tempA=getTemp(s_A);
        tempB=getTemp(s_B);
        trigger(vidobj);
        rampCurr(s_A,imin);
        rampCurr(s_B,imin);
        data=getdata(vidobj);
        [waist_x,waist_y]=guppyAnalysis_rev6(data);
        if waist_x > 1 || waist_y > 1
            keyboard
        end
        waistMatrix(indxA,indxB,1)=waist_x;
        waistMatrix(indxA,indxB,2)=waist_y;
        waistMatrix(indxA,indxB,3)=mean([waist_x waist_y]);
        tempMatrix(indxA,indxB,1)=tempA;
        tempMatrix(indxA,indxB,2)=tempB;
        indx=indx+1;
    end
end

resetCurr(s_A,s_B);
cal=waistMatrix(:,:,3);

[M,I]=min(waistMatrix(:,:,3));

uplimit=1;

while 1
    if I(uplimit)<I(uplimit+1) && I(uplimit)<3
        break
    end
    uplimit=uplimit+1;
end

M=transpose(M(1:uplimit));
xcal=transpose(currVector(1:uplimit));

%% SAVE DATA TO FILES
str=[date '\' basename '_' num2str (imin) '_' num2str (imax) '_' num2str(ppl)];
mypath = ['C:\Users\Bolton\Documents\MATLAB\Optotune\r uns\' str];
mkdir(mypath);
waistname = '_waists.csv';
tempname = '_temps.csv';

csvwrite(fullfile(mypath,['avg' waistname]),waistMatrix(:,:,3));
csvwrite(fullfile(mypath,['x' waistname]),waistMatrix(:,:,1));
csvwrite(fullfile(mypath,['A' tempname]),tempMatrix(:,:,1));
csvwrite(fullfile(mypath,['y' waistname]),waistMatrix(:,:,2));
csvwrite(fullfile(mypath,['B' tempname]),tempMatrix(:,:,2));

%% PLOTS
figure
pcolor(waistMatrix(:,:,1));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'YTick', linspace(1,ppl,10));
set(gca, 'XTick', linspace(1,ppl,10));
set(gca, 'YTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
set(gca, 'XTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
xlabel('I_B (mA)')
ylabel('I_A (mA)')
title('Guppy X Spot Diameter (mm)')

savefig(fullfile(mypath,'xspot.fig'));

figure
pcolor(waistMatrix(:,:,2));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'YTick', linspace(1,ppl,10));
set(gca, 'XTick', linspace(1,ppl,10));
set(gca, 'YTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
set(gca, 'XTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
xlabel('I_B (mA)')
ylabel('I_A (mA)')
title('Guppy Y Spot Diameter (mm)')

savefig(fullfile(mypath,'yspot.fig'));

figure
pcolor(waistMatrix(:,:,3));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'YTick', linspace(1,ppl,10));
set(gca, 'XTick', linspace(1,ppl,10));
set(gca, 'YTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
set(gca, 'XTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
xlabel('I_B (mA)')
ylabel('I_A (mA)')
% caxis([0 1.5])
title('Mean of X and Y Guppy Spot Diameters (mm)')

savefig(fullfile(mypath,'avgspot.fig'));

figure
subplot(2,2,1)
pcolor(waistMatrix(:,:,1));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'YTick', linspace(1,ppl,10));
set(gca, 'XTick', linspace(1,ppl,10));
set(gca, 'YTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
set(gca, 'XTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
xlabel('I_B (mA)')
ylabel('I_A (mA)')
% caxis([0 1.5])
title('Guppy X Spot Diameter (mm)')

subplot(2,2,2)
pcolor(waistMatrix(:,:,2));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'YTick', linspace(1,ppl,10));
set(gca, 'XTick', linspace(1,ppl,10));
set(gca, 'YTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
set(gca, 'XTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
xlabel('I_B (mA)')
ylabel('I_A (mA)')
% caxis([0 1.5])
title('Guppy Y Spot Diameter (mm)')

subplot(2,2,3)
imagesc(tempMatrix(:,:,1));
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'YTick', linspace(1,ppl,10));
set(gca, 'XTick', linspace(1,ppl,10));
set(gca, 'YTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
set(gca, 'XTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
xlabel('I_B (mA)')
ylabel('I_A (mA)')
caxis([min(min(tempMatrix(:,:,1))) max(max(tempMatrix(:,:,1)))])
title('Temperature of Lens A (\[Degree]C)')

subplot(2,2,4)
imagesc(tempMatrix(:,:,2));
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottomdele
set(gca, 'YTick', linspace(1,ppl,10));
set(gca, 'XTick', linspace(1,ppl,10));
set(gca, 'YTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
set(gca, 'XTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
xlabel('I_B (mA)')
ylabel('I_A (mA)')
caxis([min(min(tempMatrix(:,:,2))) max(max(tempMatrix(:,:,2)))])
title('Temperature of Lens B (\[Degree]C)')

savefig(fullfile(mypath,'2by2.fig'));

figure
scatter(xcal,M);
hold on
grid on
xlabel('I_B (d-mA)');
ylabel('Beam Waist (mm)');
title('Uncompensated Calibration: Average Focused Spot Size')

savefig(fullfile(mypath,'avgfocus.fig'));

keyboard
end 
function compcal=getCompCalRand(s_A,s_B,tmin,tmax,imin,imax,ppl,adrift,bdrift)
%% SETUP
clearBytes(s_A);
clearBytes(s_B);

[vidobj,~]=initializeGuppy(ppl*ppl);

basename = input('Run Name: ','s');

rampCurr(s_A,imin);
rampCurr(s_B,imin);

currVector=round(linspace(imin, imax, ppl)); % imax is given on a a scale from 0 to 4095
waistMatrix=zeros(ppl,ppl,3);
tempMatrix=zeros(ppl,ppl,2);

randindx=transpose(randperm(ppl*ppl));
refMatrix=zeros(ppl*ppl,2);

indx=1;
for ii=1:ppl
    for jj=1:ppl
        refMatrix(indx,1)=ii;
        refMatrix(indx,2)=jj;
        indx=indx+1;
    end
end

%% WARMUP
tempA=getTemp(s_A);
tempB=getTemp(s_B);

while (tempA < 26.5) || (tempB < 26.5)
    if tempA < 26.5
        rampCurr(s_A,3500);
    else
        rampCurr(s_A,2000);
    end
    if tempB < 26.5
        rampCurr(s_B,3500);
    else
        rampCurr(s_B,2000);
    end
    pause(1)
    tempA=getTemp(s_A);
    tempB=getTemp(s_B);
    disp(['Temperature A: ' num2str (tempA) ' \[Degree]C']);
    disp(['Temperature B: ' num2str (tempB) ' \[Degree]C']);
end
rampCurr(s_A,imin);
rampCurr(s_B,imin);

%% DATA COLLECTION
tic
indx=1;
for ii=1:ppl
    for jj=1:ppl
        tempA=getTemp(s_A);
        tempB=getTemp(s_B);
        disp(['Progress: ' num2str (100*(indx)/(ppl)^2) '%']);
        time=toc;
        disp(['Time: ' num2str (time) ' s']);
        eta=((ppl)^2-(ppl*(ii-1)+jj))*time/(ppl*(ii-1)+jj);
        [q, sec] = quorem(round(sym(eta)), sym(60));
        [hr, minutes]=quorem(sym(q), sym(60));
        disp(['Est:  ' num2str (int8(hr)) ':' num2str (int8(minutes)) ':' num2str(int8(sec))]);
        disp(' ');
        if (tempA >= tmax) || (tempB >=tmax)
            rampCurr(s_A,140);
            rampCurr(s_B,140);
            while (tempA >= tmax) || (tempB >=tmax)
                pause(5);
                tempA=getTemp(s_A);
                tempB=getTemp(s_B);
                disp(['Temperature A: ' num2str (tempA) ' \[Degree]C']);
                disp(['Temperature B: ' num2str (tempB) ' \[Degree]C']);
            end
        end
        indxA=refMatrix(randindx(ppl*(ii-1)+jj),1);
        indxB=refMatrix(randindx(ppl*(ii-1)+jj),2);
        currA=compCurrInterp(s_A,adrift,tmin,tmax,imin,imax,currVector(indxA));
        currB=compCurrInterp(s_B,bdrift,tmin,tmax,imin,imax,currVector(indxB));
        rampCurr(s_A,currA);
        rampCurr(s_B,currB);
        pause(0.05)
        tempA=getTemp(s_A);
        tempB=getTemp(s_B);
        trigger(vidobj);
        rampCurr(s_A,imin);
        rampCurr(s_B,imin);
        data=getdata(vidobj);
        [waist_x,waist_y]=guppyAnalysis_rev6(data);
        waistMatrix(indxA,indxB,1)=waist_x;
        waistMatrix(indxA,indxB,2)=waist_y;
        waistMatrix(indxA,indxB,3)=mean([waist_x waist_y]);
        tempMatrix(indxA,indxB,1)=tempA;
        tempMatrix(indxA,indxB,2)=tempB;
        indx=indx+1;
    end
end

resetCurr(s_A,s_B);
compcal=waistMatrix(:,:,3);

[M,I]=min(waistMatrix(:,:,3));

uplimit=1;

while 1
    if I(uplimit)<I(uplimit+1) && I(uplimit)<3
        break
    end
    uplimit=uplimit+1;
end

M=transpose(M(1:uplimit));
xcal=transpose(currVector(1:uplimit));

%% SAVE DATA TO FILES
str=[date '\' basename '_' num2str (imin) '_' num2str (imax) '_' num2str(ppl)];
mypath = ['C:\Users\Bolton\Documents\MATLAB\Optotune\r uns\' str];
mkdir(mypath);
waistname = '_waists.csv';
tempname = '_temps.csv';

csvwrite(fullfile(mypath,['avg' waistname]),waistMatrix(:,:,3));
csvwrite(fullfile(mypath,['x' waistname]),waistMatrix(:,:,1));
csvwrite(fullfile(mypath,['A' tempname]),tempMatrix(:,:,1));
csvwrite(fullfile(mypath,['y' waistname]),waistMatrix(:,:,2));
csvwrite(fullfile(mypath,['B' tempname]),tempMatrix(:,:,2));

%% PLOTS
figure
pcolor(waistMatrix(:,:,1));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'YTick', linspace(1,ppl,10));
set(gca, 'XTick', linspace(1,ppl,10));
set(gca, 'YTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
set(gca, 'XTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
xlabel('I_B (mA)')
ylabel('I_A (mA)')
% caxis([min(min(waistMatrix(:,:,1))) max(max(waistMatrix(:,:,1)))])
title('Guppy X Spot Diameter (mm)')

savefig(fullfile(mypath,'xspot.fig'));

figure
pcolor(waistMatrix(:,:,2));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'YTick', linspace(1,ppl,10));
set(gca, 'XTick', linspace(1,ppl,10));
set(gca, 'YTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
set(gca, 'XTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
xlabel('I_B (mA)')
ylabel('I_A (mA)')
% caxis([min(min(waistMatrix(:,:,1))) max(max(waistMatrix(:,:,1)))])
title('Guppy Y Spot Diameter (mm)')

savefig(fullfile(mypath,'yspot.fig'));

figure
pcolor(waistMatrix(:,:,3));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'YTick', linspace(1,ppl,10));
set(gca, 'XTick', linspace(1,ppl,10));
set(gca, 'YTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
set(gca, 'XTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
xlabel('I_B (mA)')
ylabel('I_A (mA)')
% caxis([0 1.5])
title('Mean of X and Y Guppy Spot Diameters (mm)')

savefig(fullfile(mypath,'meanspot.fig'));

figure
subplot(2,2,1)
pcolor(waistMatrix(:,:,1));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'YTick', linspace(1,ppl,10));
set(gca, 'XTick', linspace(1,ppl,10));
set(gca, 'YTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
set(gca, 'XTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
xlabel('I_B (mA)')
ylabel('I_A (mA)')
% caxis([0 1.5])
title('Guppy X Spot Diameter (mm)')

subplot(2,2,2)
pcolor(waistMatrix(:,:,2));
shading interp
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'YTick', linspace(1,ppl,10));
set(gca, 'XTick', linspace(1,ppl,10));
set(gca, 'YTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
set(gca, 'XTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
xlabel('I_B (mA)')
ylabel('I_A (mA)')
% caxis([0 1.5])
title('Guppy Y Spot Diameter (mm)')

subplot(2,2,3)
imagesc(tempMatrix(:,:,1));
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'YTick', linspace(1,ppl,10));
set(gca, 'XTick', linspace(1,ppl,10));
set(gca, 'YTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
set(gca, 'XTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
xlabel('I_B (mA)')
ylabel('I_A (mA)')
caxis([min(min(tempMatrix(:,:,1))) max(max(tempMatrix(:,:,1)))])
title('Temperature of Lens A (\[Degree]C)')

subplot(2,2,4)
imagesc(tempMatrix(:,:,2));
colorbar;
set(gca,'XDir','normal');% X origin is left
set(gca,'YDir','normal');% Y origin is bottom
set(gca, 'YTick', linspace(1,ppl,10));
set(gca, 'XTick', linspace(1,ppl,10));
set(gca, 'YTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
set(gca, 'XTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
xlabel('I_B (mA)')
ylabel('I_A (mA)')
caxis([min(min(tempMatrix(:,:,2))) max(max(tempMatrix(:,:,2)))])
title('Temperature of Lens B (\[Degree]C)')

savefig(fullfile(mypath,'2by2.fig'));

figure
scatter(xcal,M);
hold on
grid on
xlabel('I_B (d-mA)');
ylabel('Beam Waist (mm)');
title('Compensated Calibration: Average Focused Spot Size')

savefig(fullfile(mypath,'avgfocus.fig'));

keyboard
end 
  
%   III. CALIBRATION FIT FUNCTIONS: CURVE FITTING TO FOLLOW POINTS OF FOCUS
%
%     1. Here we look at the (un)compensated calibration data sets
%       collected in II. We follow the minimum line, generating a fit
%       function I_A(I_B) that follows the minima; given a requested I_B, the
%       fit function chooses I_A such that the resulting spot size is a
%       minimum.
function [focusobject, plotobject]=focusFit(picMatrix,imin,imax)
dim=size(picMatrix);
[~,I]=min(picMatrix);
x=linspace(imin,imax,dim(1));
y=imin+(imax-imin)*(I-1)/(dim(1)-1);

uplimit=1;

while 1
    if I(uplimit)<I(uplimit+1) && I(uplimit)<3
        break
    end
    uplimit=uplimit+1;
end

x=transpose(x(1:uplimit));
y=transpose(y(1:uplimit));
p=transpose(1:uplimit);
q=transpose(I(1:uplimit));

focusobject=fit(x,y,'poly3');
plotobject=fit(p,q,'poly3');

figure
imagesc(picMatrix)
hc=colorbar;
hc.Label.String='Average Spot Diameter (mm)';
hold on
set(gca,'XDir','normal');
set(gca,'YDir','normal');
set(gca, 'YTick', linspace(1,dim(2),10));
set(gca, 'XTick', linspace(1,dim(1),10));
set(gca, 'YTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
set(gca, 'XTickLabel', round(linspace(296*imin/4096,296*imax/4096,10)));
title('Mean of X and Y Spot Diameters');

h2=scatter(p,q,'black');
h3=plot(plotobject);

legend('Minima','Fit Function');
xlabel('I_B (mA)');
ylabel('I_A (mA)');
caxis([min(min(picMatrix)) max (max(picMatrix))/10]);
colorbar;
set(h3,'linewidth',2);
set(h2,'MarkerFaceColor','white');

keyboard

end
  
%   IV. CALIBRATION TESTING: UNCOMPENSATED AND COMPENSATED
%
%     1. Simple testing of (un)compensated calibration. Here we run through
%       the motions of following the minimum line without taking any data.
%       These functions are good for demonstrations.

function uncompCalTest(s_A,s_B,focusobject,imin,imax,steps,delay)
% focusobject is generated from focusFit or focusSpline

clearBytes(s_A);
clearBytes(s_B);

x=round(linspace(imin,imax,steps));
y=round(focusobject(x));

%% WARMUP
tempA=getTemp(s_A);
tempB=getTemp(s_B);

while (tempA < 26.5) || (tempB < 26.5)
    if tempA < 26.5
        rampCurr(s_A,3500);
    else
        rampCurr(s_A,2000);
    end
    if tempB < 26.5
        rampCurr(s_B,3500);
    else
        rampCurr(s_B,2000);
    end
    pause(1)
    tempA=getTemp(s_A);
    tempB=getTemp(s_B);
end
rampCurr(s_A,imin);
rampCurr(s_B,imin);

%% TEST

for ii=1:steps
    faroCurr(s_A,s_B,y(ii),x(ii),10);
    pause(delay)
end

rampCurr(s_A,0);
rampCurr(s_B,0);

end
% function compCalTest: Needs to be remade following the transition from
%       compCurrent to compCurrInterp

%     2. Same as before but now we're taking data to see how well we're
%       following the intended focsued spot sizes.
function uncompCalTestDataRand(s_A,s_B,uncompcal,uncompfocusobject,imin,imax,steps,delay)
%% SETUP
clearBytes(s_A);
clearBytes(s_B);

[vidobj,~]=initializeGuppy(steps);

basename = input('uncompCalTestDataRand Run Name: ','s');

dim=size(uncompcal);
[M,I]=min(uncompcal);
xcal=round(linspace(140,4000,dim(1)));
waistMatrix=zeros(steps,3);
tempMatrix=zeros(steps,2);

x=round(linspace(imin,imax,steps));
y=round(uncompfocusobject(x));

randindx=transpose(randperm(steps));

%% WARMUP
tempA=getTemp(s_A);
tempB=getTemp(s_B);

while (tempA < 26.5) || (tempB < 26.5)
    if tempA < 26.5
        rampCurr(s_A,3500);
    else
        rampCurr(s_A,2000);
    end
    if tempB < 26.5
        rampCurr(s_B,3500);
    else
        rampCurr(s_B,2000);
    end
    pause(1)
    tempA=getTemp(s_A);
    tempB=getTemp(s_B);
    disp(['Temperature A: ' num2str (tempA) ' \[Degree]C']);
    disp(['Temperature B: ' num2str (tempB) ' \[Degree]C']);
end
rampCurr(s_A,imin);
rampCurr(s_B,imin);

%% DATA COLLECTION
for ii=1:steps
    indx=randindx(ii);
    faroCurr(s_A,s_B,y(indx),x(indx),10);
    pause(0.05)
    tempA=getTemp(s_A);
    tempB=getTemp(s_B);
    trigger(vidobj);
    data=getdata(vidobj);
    [waist_x,waist_y]=guppyAnalysis_rev6(data);
    waistMatrix(indx,1)=waist_x;
    waistMatrix(indx,2)=waist_y;
    waistMatrix(indx,3)=mean(waistMatrix(indx,1:2));
    tempMatrix(indx,1)=tempA;
    tempMatrix(indx,2)=tempB;
    disp(['Progress: ' num2str (100*ii/steps) '%']);
    pause(delay)
end

rampCurr(s_A,0);
rampCurr(s_B,0);

uplimit=1;

while 1
    if I(uplimit)<I(uplimit+1) && I(uplimit)<3
        break
    end
    uplimit=uplimit+1;
end

M=transpose(M(1:uplimit));
xcal=transpose(xcal(1:uplimit));
calminfit=fit(xcal,M,'poly2');

%% SAVE DATA TO FILES
str=[date '\' basename '_' num2str (imin) '_' num2str (imax) '_' num2str (steps) '_' num2str(delay)];
mypath = ['C:\Users\Bolton\Documents\MATLAB\Optotune\r uns\' str];
mkdir(mypath);
csvwrite(fullfile(mypath,'atempdata.csv'),tempMatrix(:,1));
csvwrite(fullfile(mypath,'btempdata.csv'),tempMatrix(:,2));
csvwrite(fullfile(mypath,'x_waists.csv'),waistMatrix(:,1));
csvwrite(fullfile(mypath,'y_waists.csv'),waistMatrix(:,2));
csvwrite(fullfile(mypath,'avg_waists.csv'),waistMatrix(:,3));

%% PLOTS
figure
scatter(xcal,M);
hold on
grid on
scatter(x,waistMatrix(:,3));
plot(calminfit);
legend('Calibrated Focus','Cal Test: Avg Waist');
xlabel('I_B (d-mA)');
ylabel('Beam Waist (mm)');
title('Calibration Test')

savefig(fullfile(mypath,'meanfocusedspot.fig'));
keyboard
end
function compCalTestDataRand(s_A,s_B,cal,focusobject,adrift,bdrift,imin,imax,steps,delay)
%% SETUP
clearBytes(s_A);
clearBytes(s_B);

[vidobj,~]=initializeGuppy(steps);

basename = input('compCalTestDataRand Run Name: ','s');

dim=size(cal);
[M,I]=min(cal);
xcal=round(linspace(140,4000,dim(1)));
waistMatrix=zeros(steps,3);
tempMatrix=zeros(steps,2);

x=round(linspace(imin,imax,steps));
y=round(focusobject(x));

randindx=transpose(randperm(steps));

%% WARMUP
tempA=getTemp(s_A);
tempB=getTemp(s_B);

while (tempA < 26.5) || (tempB < 26.5)
    if tempA < 26.5
        rampCurr(s_A,3500);
    else
        rampCurr(s_A,2000);
    end
    if tempB < 26.5
        rampCurr(s_B,3500);
    else
        rampCurr(s_B,2000);
    end
    pause(1)
    tempA=getTemp(s_A);
    tempB=getTemp(s_B);
    disp(['Temperature A: ' num2str (tempA) ' \[Degree]C']);
    disp(['Temperature B: ' num2str (tempB) ' \[Degree]C']);
end
rampCurr(s_A,imin);
rampCurr(s_B,imin);

%% DATA COLLECTION
for ii=1:steps
    indx=randindx(ii);
    currA=compCurrInterp(s_A,adrift,25,30,140,4000,y(indx));
    currB=compCurrInterp(s_B,bdrift,25,30,140,4000,x(indx));
    faroCurr(s_A,s_B,currA,currB,10);
    pause(0.05)
    tempA=getTemp(s_A);
    tempB=getTemp(s_B);
    trigger(vidobj);
    data=getdata(vidobj);
    [waist_x,waist_y]=guppyAnalysis_rev6(data);
    waistMatrix(indx,1)=waist_x;
    waistMatrix(indx,2)=waist_y;
    waistMatrix(indx,3)=mean(waistMatrix(indx,1:2));
    tempMatrix(indx,1)=tempA;
    tempMatrix(indx,2)=tempB;
    disp(['Progress: ' num2str (100*ii/steps) '%']);
    pause(delay)
end

rampCurr(s_A,0);
rampCurr(s_B,0);

uplimit=1;

while 1
    if I(uplimit)<I(uplimit+1) && I(uplimit)<3
        break
    end
    uplimit=uplimit+1;
end

M=transpose(M(1:uplimit));
xcal=transpose(xcal(1:uplimit));
calminspline=fit(xcal,M,'poly2');

%% SAVE DATA TO FILES
str=[date '\' basename '_' num2str (imin) '_' num2str (imax) '_' num2str (steps) '_' num2str(delay)];
mypath = ['C:\Users\Bolton\Documents\MATLAB\Optotune\r uns\' str];
mkdir(mypath);
csvwrite(fullfile(mypath,'atempdata.csv'),tempMatrix(:,1));
csvwrite(fullfile(mypath,'btempdata.csv'),tempMatrix(:,2));
csvwrite(fullfile(mypath,'x_waists.csv'),waistMatrix(:,1));
csvwrite(fullfile(mypath,'y_waists.csv'),waistMatrix(:,2));
csvwrite(fullfile(mypath,'avg_waists.csv'),waistMatrix(:,3));

%% PLOTS
figure
scatter(xcal,M);
hold on
grid on
scatter(x,waistMatrix(:,3));
plot(calminspline);
legend('Calibrated Focus','Cal Test: Avg Waist');
xlabel('I_B (d-mA)');
ylabel('Beam Waist (mm)');
title('Calibration Test')

savefig(fullfile(mypath,'meanfocusedspot.fig'));

keyboard
end
