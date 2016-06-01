function [vidobj,src]=initializeGuppy(datapoints,exposure)
vidobj = videoinput('gentl', 1);
src = getselectedsource(vidobj);

triggerconfig(vidobj,'manual');
vidobj.FramesPerTrigger=1;
vidobj.TriggerRepeat=datapoints-1;
% vidobj.ROIPosition = [250 150 148 148];
src.ExposureTime=exposure;

start(vidobj);
end

