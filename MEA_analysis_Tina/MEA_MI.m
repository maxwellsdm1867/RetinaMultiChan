close all
clear all

%%%%%%%%%%%%%%%%%%%%% call data %%%%%%%%%%%%%%%%%%%%%%%%
% filename = 'E:\Tina\data_MEA\G15-3_20120309_20120402_CM_0002.mcd'
[filename,directory] = uigetfile('*.mcd','Pick the AXON file'); 
file = strcat(directory,filename)
%%% determine the channel length
data1 = mreadMCD(file, MEA_layout(1));
data2 = mreadMCD(file, MEA_layout(2)); ;
h = mutualinfo(data1,data2) ;
% %%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%
% n_rows=length(data);
% SamplingRate = 20000.;
% t_data = [1:n_rows] / SamplingRate ;
% DataTime = n_rows / SamplingRate ;
% 
% %%%%%%%%%%%%%%%%%%%% read mcd data channel by channel  %%%%%%%%%%%%%%%%%
% for channel_index = 1 : 60
%     channel_index
%     data1 = mreadMCD(file, MEA_layout(channel_index));   
% end
% for channel_index = 1 : 60
%     channel_index
%     data2 = mreadMCD(file, MEA_layout(channel_index));   
% end

