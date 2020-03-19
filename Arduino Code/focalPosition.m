%This is just some preliminary code to work with getting power
%and intensity reading from the photodiodes with the Arduino
%I will also try to get the focal position information with this
%code before it is used to feedback to setting the focal power of
%the lenses.
function focalPosition
COM='COM8';
timeLimit = 100;
numPoints = 1000;
dataCounter = 1;

try
    delete(instrfind('Port',COM));
    s=establishConnection(COM);
    
    while(s.BytesAvailable>0)
        fread(s,s.BytesAvailable);
    end
    pause(0.5);
catch
    delete((instrfind('Port',COM)));
    error('Could not connect to Arduinos.');     
end



mygui=figure('Name' ,'Seed Monitor', ...
    'Position', [100 400 1000 600], ...
    'CloseRequestFcn', @closeFcn, ...
    'Visible', 'Off');
set(mygui,'Resize','Off');
set(mygui,'MenuBar','none');
set(mygui,'NumberTitle','Off');
set(mygui,'Color','w');
mygui.Visible = 'On';


voltageaxes = subplot(2,1,1);
title('Photodiode Voltage vs Time')
xlabel('Time')
ylabel('Voltage')
set(voltageaxes, 'XMinorGrid', 'On', 'YMinorGrid', 'On', ...
    'XMinorTick', 'On', 'YMinorTick', 'On', 'Box', 'On');
set(voltageaxes,'FontSize',8);

hold(voltageaxes,'on');

powers = scatter(voltageaxes,0,0);
powers.Marker='o';
powers.MarkerFaceColor='r';
powers.MarkerEdgeColor='r';
powers.SizeData = 50;

intensities = scatter(voltageaxes,0,0);
intensities.Marker='o';
intensities.MarkerFaceColor='k';
intensities.MarkerEdgeColor='k';
intensities.SizeData = 50;

focalaxes = subplot(2,1,2);
title('Calculated Focal Position vs. Time (at a constant offset)')
xlabel('Time')
ylabel('Focal Position')
set(focalaxes, 'XMinorGrid', 'On', 'YMinorGrid', 'On', ...
    'XMinorTick', 'On', 'YMinorTick', 'On', 'Box', 'On');
set(focalaxes,'FontSize',8);

hold(focalaxes,'on');

foci = scatter(focalaxes,0,0);
foci.Marker='o';
foci.MarkerFaceColor='b';
foci.MarkerEdgeColor='b';
foci.SizeData = 50;

GOAT = legend(voltageaxes,{'Total Power','Central Intensity'},'Location','northeast')

tic;

while (dataCounter < numPoints)
    [power,intensity] = readValue();
    powers.XData(dataCounter) = toc;
    powers.YData(dataCounter) = power;
    intensities.XData(dataCounter) = toc;
    intensities.YData(dataCounter) = intensity;
    foci.XData(dataCounter) = toc;
    focus = findFocus(power,intensity);
    foci.YData(dataCounter) = focus;
    dataCounter = dataCounter + 1;
    pause(0.05);
end


fclose(s);

function closeFcn(~,~)
    disp('Closing GUI...');
    disp('Closing the serial connection.');
    fclose(s);
    disp(['Deleting serial connections on PORT ' COM]);
    delete(instrfind('Name', ['Serial-' COM]));
    disp('Deleting the figure...');
    delete(mygui);    
end





% function value=readValue
%     fwrite(s,10);               % write a line feed
%     value=fscanf(s);            % read it
%     value(end)=[];              % Get rid of lf charcter
%     value=str2double(value);     % Convert to number
% end

function [power,intensity] = readValue
    fwrite(s,13);
    power = fscanf(s);
    power(end) = '';
    power = str2double(power);
    intensity = fscanf(s);
    intensity(end) = '';
    intensity = str2double(intensity);
%     disp(output)
end

function position = findFocus(power,intensity)
    %power readings: 4.3 into PBS
    %0.7 into power photodiode
    %NOT MEASURED DIRECTLY but stands to reason ~3.6 in intensity
    %so to make sure that things are on equal footing in the beam
    %we should multiply the power reading by (4.3/0.7)
    %and multiply the intensity reading by (4.3/3.6)
    pratio = 4.3/0.6;
    iratio = 4.3/3.6;
    waist = (1.5E-4/2); %meters
    lambda = 1.064E-6; %meters
    rayleighrange = (pi*waist^2)/(lambda); %meters
    irisradius = 7.5E-4; %it's all in meters!!
    irisarea = pi*(irisradius^2);
    power = power*pratio;
    intensity = intensity*iratio;
    intensity = intensity/irisarea;
    ratio = intensity/power;
    position = rayleighrange*sqrt(2/(pi*ratio*(waist^2))-1)-0.1;
    %now that these values are on equal footing we can comapre
    %them using the typical gaussian beam equations
    
end

end



function s=establishConnection(COM)
    port=COM;
    baud=9600;
    databits=8;
    stopbits=1;
    parity='none';
    flowcontrol='none';
    timeout = 0.5;
    s=serial(port,'BaudRate',baud,...
        'Parity',parity,'StopBits',...
        stopbits,'DataBits',databits,'FlowControl',flowcontrol,'Timeout',...
        timeout);
    set(s,'terminator','CR');

    fopen(s);
    disp('Connection has been establishied.');
end

