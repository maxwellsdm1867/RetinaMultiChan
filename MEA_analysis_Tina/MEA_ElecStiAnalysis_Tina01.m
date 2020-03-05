clear all
close all

%%%%%%%%%%%%%%% simple code for stella %%%%%%%%%%%%%%%%%%
% using butter filter to remove shift
file = 'E:\MC_Rack Data\20140512\20140418(MEA(PtA2-9))_20140512001.mcd'
% file = strcat(directory,filename)
% filenameold = filename ;
%%% determine the channel length
AllDataInfo =datastrm(file) ;
%if ~exist('t1','var')
    t1=getfield(AllDataInfo,'sweepStartTime');
%end
%if ~exist('t2','var')
    t2=getfield(AllDataInfo,'sweepStopTime') ;
%end
StartStopVector = [t1 t2];
AllData = nextdata(AllDataInfo,'startend',StartStopVector,'originorder','on' ,'streamname','Electrode Raw Data');
rawdata = AllData.data(1,:);
SamplingRate = getfield(AllDataInfo,'MillisamplesPerSecond2')/1000 ; % unit is Hz
%n_orgData = length(data) ;
%data = data(t1 * SamplingRate + 1 : t2 * SamplingRate);

%%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%
n_rows=length(rawdata);
DataTime = n_rows / SamplingRate ;

%%%%%%%%%%%%%%%%%%%% read mcd data channel by channel  %%%%%%%%%%%%%%%%%

for channel_index = 1 : 60
    tic
    channel_index
    rawdata = AllData.data(MEA_layout_Tina02(channel_index),:); 
    [vdata, tvals] = ad2muvolt(AllDataInfo, rawdata,  getfield(AllDataInfo,'StreamNames'));
%     data = vdata' ;
%     dataBaseLine = smooth(data,1000) ;
%     data = data - dataBaseLine ;
    [b,a] = butter(2,200/10000,'high'); % set butter filter (order, threshold/(sampling rate/2),high pass)
    FilterData = filter(b,a,vdata);
    %%% peak detecttion
    mean_data = mean(FilterData) ;
    std_data = std(FilterData) ;
    [pks,locs] = findpeaks(-FilterData,'MINPEAKHEIGHT', std_data*4);
%     p = peakfinder(data, 0, -1*std_data*5,-1);
    plot(-FilterData), hold on, plot(locs,std_data*5,'o'),hold off 
    SpikeTime{channel_index} = locs/SamplingRate ;
    toc
end


