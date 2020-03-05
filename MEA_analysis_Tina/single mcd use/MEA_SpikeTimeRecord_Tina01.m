close all
clear all

%%%%%%%%%%%%%%%%%%%%% call data %%%%%%%%%%%%%%%%%%%%%%%%
% filename = 'E:\Tina\data_MEA\G15-3_20120309_20120402_CM_0002.mcd'
[filename,directory] = uigetfile('*.mcd','Pick the AXON file'); 
file = strcat(directory,filename)
%%% determine the channel length
data = mreadMCD(file, MEA_layout(1));

%%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%
n_rows=length(data);
SamplingRate = 20000.;
t_data = [1:n_rows] / SamplingRate ;
DataTime = n_rows / SamplingRate ;

%%%%%%%%%%%%%%%%%%%% read mcd data channel by channel  %%%%%%%%%%%%%%%%%
for channel_index = 1 : 60
    tic
    channel_index
    data = mreadMCD(file, MEA_layout(channel_index));   
    
    %%% peak detecttion
    mean_data = mean(data) ;
    std_data = std(data) ;
    p = peakfinder(data, -1, -1*std_data*5,-1);
%     plot(t_data(1:1000000), data(1:1000000) , p(1:100)/SamplingRate ,std_data*-5 , 'o') 

    SpikeTime{channel_index} = p/SamplingRate ;
    toc
end

filename = ['MEA_spike_' filename(1:length(filename)-4) '.mat']
file = strcat(directory,filename) 
save(file,'SpikeTime')
'finish MEA_SpikeTimeRecord_Tina01'