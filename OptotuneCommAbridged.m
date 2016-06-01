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

keyboard

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
%This function clears any serial communication bytes that are sitting in
%channel s
while s.BytesAvailable>0
    fread(s,1);
end
end

function s=establishConnection(portno,portname)
%This function establishes a serial connection between MATLAB and the
%desired com port
port=['COM' num2str(portno)];
baud=115200;
databits=8;
stopbits=1;
parity='none';
flowcontrol='none';

display(sprintf('Opening %s...',port));
s = serial(port,'BaudRate',baud,'Parity',parity,'StopBits',stopbits,'DataBits'
	,databits,'FlowControl',flowcontrol,'Name',portname);

set(s,'terminator','CR');
fopen(s);

disp(s);
disp('Serial connection established.');
disp(' ');

end

function AO=generateAO
%This function uses the MCC USB-1408FS-Plus and creates an analog output
%object (AO). It then adds the two analog output channels (0 and 1) and
%renames them a_A and a_B (analog output for driver A and B, respectively).
%The second parameter on analogoutput() might need to be changed depending
%on whether or not there are other MCC devices on the computer.

AO=analogoutput('mcc',1);
addchannel(AO,0:1,{'a_A','a_B'});
AO.SampleRate=25000;
disp('USB-1408FS-Plus: AO device ready.')
disp(' ');
end

% GENERAL READ/WRITE COMMANDS: Talk to the lens drivers and DAQ.
function handshake(s)
%Is the Lens Driver available?
clearBytes(s)
fprintf(s,'%s\r','Start');
reply=fread(s,7);
disp(['Lens ' get(s,'Name') ' ' char(transpose(reply(1:5)))]);
end

function getID(s)
%This gets the ID of the Lens Driver connected to port s. This also
%verifies that the driver is available.
clearBytes(s)
mystr='IRaaaaaaaa';%5993
bin=uint8(mystr);
high=uint8(hex2dec('59'));
low=uint8(hex2dec('93'));

msg=[bin low high];

fwrite(s,msg);
reply=fread(s,14);

disp(['Board ID ' get(s,'Name') ': ' char(transpose(reply(3:10)))]);
end

function getFirmware(s)
%Is the driver up to date?

clearBytes(s)
% mystr='Va';%48FE
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

disp(['Firmware Version ' get(s,'Name') ': ' a '.' b '.' c '.00' d]);
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

disp(['Status ' get(s,'Name') ':'])

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
%What is the temperature of the lens connected to the lens driver on port s?
clearBytes(s)
fwrite(s,[84    67    65   176   208]);
reply=fread(s,9);
temp=bin2dec([dec2bin(reply(4),8) dec2bin(reply(5),8)])/16;
% disp(['Temperature ' get(s,'Name') ': ' num2str(temp) ' degC']);
% disp(' ');
end

function alog=getAnalog(s) %#ok<*DEFNU>
%This checks the ADC on the lens driver connected to s. The result is given
%as a 10 bit number 4*(0-1023), mapped from 0-5 V.
clearBytes(s)
%mystr=double('GAA');
%msg=appendCRC16(mystr);

fwrite(s,[71 65 65 64 117]);
reply=fread(s,7);
alog=4*bin2dec([dec2bin(reply(4),8),dec2bin(reply(5),8)]);
end

function volt=getVolt(s)
%Same as getAnalog(s) but now we're converting the measurement to volts.
clearBytes(s)
%mystr=double('GAA');
%msg=appendCRC16(mystr);

fwrite(s,[71 65 65 64 117]);
reply=fread(s,7);
out=bin2dec([dec2bin(reply(4),8),dec2bin(reply(5),8)]);
volt=5*1.0170*(out/1024);
end

function setAnalog(AO,atarget_new,btarget_new)
%This immediately pushes the 12 bit values atarget and btarget to
%channels a_A and a_B on the analog output device AO. Cal is:
%a_A: 5V in -> 991/1023=4.8436 V out
%a_B: 5V in -> 993/1023=4.8534 V out

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

asteps=5*transpose(linspace(atarget_old,atarget_new,100))/4095;
bsteps=5*transpose(linspace(btarget_old,btarget_new,100))/4095;
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
        disp(['Mode ' get(s,'Name') ': ' char(modeArray(ii,:))]);
    end
end

disp(' ')
end

function modeDC(s)
%This function switches the lens to DC mode from any other mode. It seems
%to work as modeSerial is supposed to work when switching from Analog mode.

% mystr=double('MwDA');
% msg=appendCRC16(mystr);
fwrite(s,[77   119    68    65    84    70]);
clearBytes(s)
end

function modeAnalog(s)
%This switches the lens to Analog mode where the current is controlled by
%the 10 bit ADC on the lens driver.
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
%In theory this brings us back to serial control of the currents in the
%lenses, or at a minimum disables analog control. This might be useless as
%modeDC does a better job of serving this purpose.
fwrite(s,[76    83   116   253]);
clearBytes(s)
end

% PROPERTY READ/CHANGE COMMANDS: Manage operating mode parameters.
function high=getHigh(s)
%This gets the upper bound for current swings in the three periodic lens modes.
clearBytes(s);
% mystr=double('PrUAaaaa');
% msg=appendCRC16(mystr);
fwrite(s,[80   114    85    65    97    97    97    97   193   107]);
reply=fread(s,9);
high=bin2dec([dec2bin(reply(4),8) dec2bin(reply(5),8)]);
end

function low=getLow(s)
%This gets the lower bound for current swings in the three periodic lens modes.
clearBytes(s);
% mystr=double('PrLAaaaa');
% msg=appendCRC16(mystr);
fwrite(s,[80   114    76    65    97    97    97    97   195    98]);
reply=fread(s,9);
low=bin2dec([dec2bin(reply(4),8) dec2bin(reply(5),8)]);
end

function freq=getFreq(s)
%This gets the frequency for the three periodic lens modes.
clearBytes(s);
% mystr=double('PrFAaaaa');
% msg=appendCRC16(mystr);
fwrite(s,[80   114    70    65    97    97    97    97   195   200]);
reply=fread(s,11);
freq=bin2dec([dec2bin(reply(6),8) dec2bin(reply(7),8)])/1000;
end

function setHigh(s,high)
%This sets the upper bound for current swings in the three periodic lens modes.
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
%This sets the lower bound for current swings in the three periodic lens modes.

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
%See getFreq. Change set -> get. ALSO: Frequency has a lower bound of 0.06
%Hz, an upper bound of 50 Hz. Setting any lower frequency will default to 1
%Hz, and the upper bound was imposed to protect the lens.

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
%This requests the maximum current (imposed by the hardware-- DO NOT CHANGE
%THIS), and the lower and upper current bounds (imposed by the software--
%currently maximized; really no reason to make it any smaller). 

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
max=bin2dec([dec2bin(maxreply(4),8) dec2bin(maxreply(5),8)])/100;

clearBytes(s)
fwrite(s,[67   114    76    65    97    97   152     4]);
lowreply=fread(s,9);
lower=max*bin2dec([dec2bin(lowreply(4),8) dec2bin(lowreply(5),8)])/4095;

clearBytes(s)
fwrite(s,[67   114    85    65    97    97   159    88]);
uppreply=fread(s,9);
upper=max*bin2dec([dec2bin(uppreply(4),8) dec2bin(uppreply(5),8)])/4095;
end

% CURRENT READ/WRITE COMMANDS: Change the current.
function curr=getCurr(s)
%Gets the current (reported as a 12 bit number 0-4096).

clearBytes(s)
fwrite(s,[65   114     0     0   180    39]);
reply=fread(s,7);
curr=bin2dec([dec2bin(reply(2),8) dec2bin(reply(3),8)]);
end

function setCurr(s,curr)
%Sets the current. Curr is a 12 bit number, NOT in mA.

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
%Ramps the current in ONE of the lenses between the present current and the
%target current. Not very smooth, but it's better than simply setting the
%current to the target. This process is done via SERIAL and is thus not
%very efficient. Use the analog control method for smoother, albeit more
%roundabout, results.
curr=getCurr(s);

currentVector=round(linspace(curr,target,10));

for ii=1:length(currentVector)
    setCurr(s,currentVector(ii));
end
end

function faroCurr(s_A,s_B,atarget,btarget,res)
%Same thing as rampCurr but this is designed to ramp the currents in both
%lenses simultaneously (Faro shuffle: Google it).

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
%When we're done and we need to semi-smoothly stop the current in the lenses.

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