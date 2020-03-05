function [mcd_data, file, channels] = mreadMCD(file, channels, t1, t2)
% NOTE: This code is not modified for the arbitrary size of file, it only
% works on the ple of the 2000000 data points(10(sec.) recording)

%The code is to read the part of the raw data of the MC_Rack recording by
%#input the beginning time and end time
%input:
% file:The location of the file to be read
% channels:Which channels to be read; you can select single one, or multi-
%          select different channels by array form [ch01, ch02, ch05,...]
%          Note that the channel indice correspond to the layout of the MEA
%          array such that you need to know the corresponding labels of the 
%          channels
% t1:read the file from t1, in (sec.) 
% t2:Read the file to t2, in (sec.)
% #output:
% mcd_data:The result in an array of [m, n]; m refers to the number of 
%          channels; n refers to the points in the selected time interval
%          The values refer to the recorded potential in microvolt
%file:The location of the file which has read
%channels:The channels which are selected

% disp('reading...');
dll = 'E:\Tina\Dropbox\MatlabCode\MEA\nsMCDLibrary64_3.4\Matlab\Matlab-Import-Filter\Matlab_Interface\nsMCDLibrary64.dll';
% define the sampling rate
Fs = 20000; 

tic;
if(~(exist('file')))
    [filename, pathname] = uigetfile('.mcd','Choose the MEA file');
    file = [pathname,filename];
end

[ns_result] = ns_SetLibrary(dll);
if(ns_result ~= 0)
    disp('DLL is not found!');
    return
end

[ns_result, hfile] = ns_OpenFile(file);
if(ns_result ~= 0)
    disp('file can not be open!');
    return
end

[ns_result, fileinfo] = ns_GetFileInfo(hfile);
if(ns_result ~= 0)
    disp('file information can not be accessed!');
    return
end

[ns_result, entityinfo] = ns_GetEntityInfo(hfile, [1:fileinfo.EntityCount]);
analog_address = find([entityinfo.EntityType] == 2);
analog_num = length(analog_address);
if(analog_num == 0)
    disp('No analog data!');
    return
end
if(analog_num ~= 60)
    disp('The channels are not all active!');
end

for analog_entity_index=1:analog_num
    analog_entity_label_string = strread(entityinfo(analog_address(analog_entity_index)).EntityLabel, '%s');
    analog_label(analog_entity_index) = str2num(analog_entity_label_string{3});
    layout_label(analog_entity_index) = str2num(analog_entity_label_string{4});
end

if(~(exist('channels')))
    layout_label
    disp(['Enter the channels of interest    ';'(Entering nothing is all channels)']);
    channels = input(':');
    if(isempty(channels))
            channels = layout_label;
    end
end

channel_num = length(channels);
% mcd_data = zeros(channel_num, (t2-t1)*Fs);

for entity_index=1:channel_num
    layout_label_address=find(layout_label == channels(entity_index));
    points = entityinfo(analog_address(layout_label_address)).ItemCount;
    if(~(exist('t1')))
        t1 = 1; 
    end
    if(~(exist('t2')))
        t2 = points/20000;
    end
%     t1_points = t1*20000+1;
%     t2_points = t2*20000; 
    [ns_result, analoginfo] = ns_GetAnalogInfo(hfile, analog_address(layout_label_address));
%     smooth_N = ceil(points/2000000);    %100(sec.) to smooth
    [ns_result, continuouscount, data] = ns_GetAnalogData(hfile, analog_address(layout_label_address), 1, points);
%     [ns_result, continuouscount, data] = ns_GetAnalogData(hfile, analog_address(layout_label_address),(entity_index-1)*point-1 , points);
    mcd_data = data ;
%     for smooth_index = 1:smooth_N
%     [ns_result, continuouscount, data] = ns_GetAnalogData(hfile, analog_address(layout_label_address), (points/smooth_N)*(smooth_index-1)+1, points/smooth_N);
%     sdata(1,:) = data;
%     sdata = sdata-smooth(sdata,100)';
%     mcd_data(entity_index, ((points/smooth_N)*(smooth_index-1)+1):points/smooth_N*smooth_index) = sdata(1,:)'*1000000;%just make the voltage unit to microvolt
%     end
end
% disp('read OK');

    