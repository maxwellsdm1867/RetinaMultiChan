% 2014 Apr.  by Tina
% read *.mcd data, record sipke time by mc_rack support
% 'C:\Program Files (x86)\Multi Channel Systems\MC_Rack\MCStreamSupport\Matlab\meatools\mcintfac'
% 
% SpikeTime,  1 * 60 cell, each cell is spiking time of a channel
% MEABurst, 1*60 cell, each cell is n cells with a series of bursts
% MEABurstTime, 1*60 cell, each cell is spiking time in bursts of a channel

SamplingRate = 20000.;

%%%%%%%%%%%%%%%%%%%%%%%%%% call data %%%%%%%%%%%%%%%%%%%%%%%%
% filename = 'E:\To Jastin\Mg\F7-2\100.11.18(MEA(F7-2))_20111205002.mcd'
%[filename,directory] = uigetfile('*.mcd','Pick the AXON file'); 
file = strcat(directory,filename)
filenameold = filename ;
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
    data = vdata' ;
    dataBaseLine = smooth(data,1000) ;
    data = data - dataBaseLine ;
    %%% peak detecttion
    mean_data = mean(data) ;
    std_data = std(data) ;
    [pks,locs] = findpeaks(-data,'MINPEAKHEIGHT', std_data*4);
%     p = peakfinder(data, 0, -1*std_data*5,-1);
%     plot(-data), hold on, plot(locs,std_data*5,'o'),hold off 
    SpikeTime{channel_index} = locs/SamplingRate ;
    toc
end

filename = ['MEA_spike_' filename(1:length(filename)-4) '.mat']
file = strcat(directory,filename) 
save(file,'SpikeTime')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% binning data with 5ms interval, use function 'BinningData_Tina01'
% DataTime = input('how many seconds of this data?  ')
t_line = [0:0.005:DataTime] ;
[BinningSpike] = BinningData_Tina01(SpikeTime,DataTime) ;
imagesc(BinningSpike)

name = ['BinningSpike_' filename(10:length(filename)-4) '.mat']
file = strcat(directory,name) 
save(file,'BinningSpike')
name = ['BinningSpike_' filename(10:length(filename)-4) '.jpg']
file = strcat(directory,name)
h = mat2gray(BinningSpike) ;
imagesc(t_line+t1,[1:60],~h)
colormap('gray')
title(name)
xlabel('time / s', 'Fontsize',16)
ylabel('channel ID', 'Fontsize',16)
% F1=getframe(gcf) ; get the image from figure 'getframe(gcf)' or 'getframe'
% imwrite(F1.cdata,file,'jpg')
saveas(gcf,file,'jpg'); % get the image from figure 'gca' or 'gcf'
'finish spike binning'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% burst detect , use function 'BurstDetect_Tina01', binning spike data which are inside a burst 
for i = 1:60
    set_SpikeInterval = 0.5 ;
    set_SpikeInBurst = 5 ;
    aa = SpikeTime{i} ;
    aa = aa' ;
    [burst , BurstTime] = BurstDetect_Tina01(aa,set_SpikeInterval,set_SpikeInBurst) ;
    MEABurst{i} = burst ;
    MEABurstTime{i} = BurstTime ;
end
name = ['MEA_Burst_' filename(10:length(filename)-4) '.mat']
file = strcat(directory,name)
save(file,'MEABurstTime','MEABurst')
'finish burst detect'

[BinningBurst] = BinningData_Tina01(MEABurstTime,DataTime) ;
imagesc(BinningBurst)
name = ['BinningBurst_' filename(10:length(filename)-4) '.mat']
file = strcat(directory,name)
save(file,'BinningBurst')

name = ['BinningBurst_' filename(10:length(filename)-4) '.jpg']
file = strcat(directory,name)
h = mat2gray(BinningBurst) ;
imagesc([0:0.005:DataTime],[1:60],~h)
title(name)
xlabel('time / s', 'Fontsize',16)
ylabel('channel ID', 'Fontsize',16)
colormap('gray')
% F1=getframe ;
% saveas(F1.cdata,file,'jpg')
saveas(gcf,file,'jpg');
'finish burst binning'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% histogram analysis,  sum spikes along the channel axis ,seperate burst by BurstDetect_Tina01
name = ['BurstHist_' filename(10:length(filename)-4)]
file = [strcat(directory,name) '.tif' ] 

SumBurst = sum(BinningSpike) ;
SumBurstCutInd = find(SumBurst < 10) ;
SumBurstCut = SumBurst ;
SumBurstCut(SumBurstCutInd) = 0 ;
SumBurstTime = t_line(SumBurstCut > 0.5) ;
set_SpikeInterval = 0.1 ;
set_SpikeInBurst = 10 ;

[HistBurst , HistBurstTime] = BurstDetect_Tina01(SumBurstTime,set_SpikeInterval,set_SpikeInBurst) ;
% plot histogram
SepBurst(1).his = [] ;
SepBurst(1).time = [] ;
for i = 2 :length(HistBurst)-1
    b = HistBurst{i} ;
    bb = find(t_line == b(1)) ;
    bl = length(b) ;
    bbb = find (t_line == b(bl)) ;
    if abs(b(1) - b(length(b))) > 0.9
        SepBurst(i).his = SumBurst(bb-20:bbb+20) ;
        SepBurst(i).time = t_line(bb-20:bbb+20) ;
    else
        SepBurst(i).his = SumBurst(bb-20:bb+180) ;
        SepBurst(i).time = t_line(bb-20:bb+180) ;
    end
    bar(SepBurst(i).time,SepBurst(i).his) ;
    title(name)
    xlabel('time / s', 'Fontsize',16)
    ylabel('count', 'Fontsize',16)
    set(gcf,'color',[1 1 1]) % set gackground as white
    F1=getframe(gcf) ; 
    imwrite(F1.cdata,file,'writemode','append')
end
name = ['Histo' filename(10:length(filename)-4) '.mat']
if ~exist('Histo', 'dir')
    mkdir(directory,'Histo') 
end;
DirectorySub = [directory 'Histo\'] ;
file = [strcat(DirectorySub,name) ] 
save(file,'SepBurst')

'finish burst histgram'
