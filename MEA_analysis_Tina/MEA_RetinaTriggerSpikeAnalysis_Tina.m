%%%%%  sweeps data spike detection by Tina @ May 2014


% file = strcat(directory,filename) ; 
clear all
close all

%%%%%%%%%%%%%%%%%%%%%% call all data in %%%%%%%%%%%%%%%%%%%%%%%%%%

file = 'C:\Users\pagepage\Desktop\retina\140718\0718_with_moved_bar.mcd';

AllDataInfo =datastrm(file); % open data
n_channel = getfield(AllDataInfo,'NChannels2'); % # of channels
n_channel = n_channel(2);% # of channels

%%%%%%%%%%%%%%%%%%%% call sweeps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_sweeps = length( getfield(AllDataInfo,'sweepStartTime'))-1 ;
% figure(3)
for sweepindex = 1:n_sweeps; % choose the sweepindex
sweepindex
SamplingRate = getfield(AllDataInfo,'MillisamplesPerSecond2')/1000; % rescale sample rate to Hz
SamplingRate =  SamplingRate(2);

sweeptimes = getfield(AllDataInfo,'sweepStartTime'); %the start time of sweep
timeWindow = getfield(AllDataInfo,'TimeWindow'); % the duration of swep
startend = [sweeptimes(sweepindex) sweeptimes(sweepindex)+timeWindow.Time2]; % startend time of signal
AllData = nextdata(AllDataInfo,'startend',startend,'originorder','on' ,'streamname','Electrode Raw Data 1'); % call data, the 60 channels
rawdata = AllData.data ;

[vdata, tvals] = ad2muvolt(AllDataInfo, rawdata', 'Electrode Raw Data 1'); % rescale vdata(uV) , tvals(us)
tvals = tvals/1E6; % rescale unit to s
% plot(tvals,vdata);
vdata2D = reshape(vdata,n_channel,length(rawdata)/n_channel) ; % separate 60 channels
tvals2D = reshape(tvals,length(rawdata)/n_channel,n_channel) ;
Time = tvals2D(:,1)';
BinningInterval = 5/1000 ; % (s) binning duration = 5ms
BinningTime = [0:BinningInterval:Time(end)] ;% # 0~600s with 5ms interval

%%%%%%%%%%%%%%%%%%%%%%%% process channel by channel %%%%%%%%%%%%%%%%%%%%%
for channel_index = 1:n_channel;
%     channel_index
    
    data = vdata2D(channel_index,:)';
    [b,a] = butter(2,200/10000,'high'); % set butter filter
    FilterData = filter(b,a,data);
%     figure(2),subplot(211),plot(Time,data),subplot(212),plot(Time,FilterData)
    std_data = std(FilterData);
    [pks,locs] = findpeaks(FilterData,'MINPEAKHEIGHT', std_data*5);
%     plot(-FilterData), hold on, plot(locs,std_data*5,'o'),hold off 
    SpikeTime{channel_index} = locs/SamplingRate ;
end

% SpikeTime = SpikeTime(~cellfun('isempty',SpikeTime)); % remove  empty cell; remove the channel without spikes 

BinningSpike{sweepindex} = BinningData_Tina01(SpikeTime,Time(end)) ;% all info is here!!!!
subplot(n_sweeps+1,1,sweepindex);
plot(BinningTime,sum(BinningSpike{sweepindex}));

end

% sweeptimes = getfield(AllDataInfo,'sweepStartTime'); %the start time of sweep
% timeWindow = getfield(AllDataInfo,'TimeWindow'); % the duration of swep
% sweepindex = 1;
startend = [sweeptimes(sweepindex) sweeptimes(sweepindex)+timeWindow.Time2];
TriggerAllData = nextdata(AllDataInfo,'startend',startend,'originorder','on' ,'streamname','Analog Raw Data'); % call trigger data
TriggerData = TriggerAllData.data(2,:);
TriggerData2 = [downsample(TriggerData,100) TriggerData(end)];
TriggerTime2 = [downsample(Time,100) Time(end)];
[TriggerData3, tvals] = ad2muvolt(AllDataInfo, TriggerData2', 'Electrode Raw Data 1');

subplot(n_sweeps+1,1,sweepindex+1);
plot(TriggerTime2, TriggerData3,'r')

samexaxis('abc','xmt','on','ytac','join','yld',1) % function need ''parseArgs'' function, from matlab center
axes;
h = title('different sweeps (sum 60 channels of a sweep)');
set(gca,'Visible','off');
set(h,'Visible','on');

figure(4)
SumSweeps = sum(cat(n_sweeps,BinningSpike{:}),n_sweeps);
imagesc(BinningTime,[1:60],SumSweeps);

title('different channels (sum sweeps of a channel)');
ylabel('Channel ID');

%%%%%%%%%%%%%%%%%%%%%% latency calculaiton %%%%%%%%%%%%%%%%%%%%%%%%
DiffTrigger = diff(TriggerData3);
EndTriggerInd =find(DiffTrigger < -1000);
LastEndTriggerInd = EndTriggerInd(end)+1;
StartTriggerInd =find(DiffTrigger > 1000);
LastEndTriggerInd = StartTriggerInd(end)+1;
LastEndTriggerTime = BinningTime(LastEndTriggerInd) ;

AllSpikes = sum(SumSweeps) ;
figure
plot(BinningTime, AllSpikes)
hold on 
plot(LastEndTriggerTime,max(AllSpikes),'vr')
% plot(BinningTime(StartTriggerInd+1),max(TriggerData3)/10,'vg')
hold off
% 
% std_data = std(AllSpikes) ;
% [pks,locs] = findpeaks(AllSpikes,'MINPEAKHEIGHT', std_data);
% LastBurstTime = BinningTime(locs(end));
% Latency = LastEndTriggerTime - LastBurstTime ; 
% 
% SumBurstCutInd = find(AllSpikes < 20) ;
% SumBurstCut = AllSpikes ;
% SumBurstCut(SumBurstCutInd) = 0 ;
% SumBurstTime = BinningTime(SumBurstCut > 0.5) ;
% set_SpikeInterval = 0.1 ;
% set_SpikeInBurst = 10 ;
% 
% [HistBurst , HistBurstTime] = BurstDetect_Tina01(SumBurstTime,set_SpikeInterval,set_SpikeInBurst) ;

for i = 1 : n_channel
    dump1 = SumSweeps(i,:);
    std_data = std(dump1) ;
    [pks,locs] = findpeaks(dump1,'MINPEAKHEIGHT', std_data);
    Dump2 =  BinningTime(locs) - LastBurstTime ;
    [MinLocs,MinLocsInd] = min(abs(Dump2));
    ChannelLastBurst(i) = BinningTime(locs(MinLocsInd));
    ChannelLatency(i) = ChannelLastBurst(i) - LastEndTriggerTime ; 
end

ChannelLatency2 = ChannelLatency;
ChannelLatency2(abs(ChannelLastBurst-LastBurstTime)> 0.05)= [] ;
StdChannel = std(ChannelLatency2);

'end MEA_RetinaTriggerSpikeAnalysis_Tina'