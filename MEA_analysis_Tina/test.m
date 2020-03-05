%%% funciton same as MEA_SpikeTimeAnalysis_Tina
%%% use in MEA_SpikeTimeAnalysis_Folderuse program





%%%%%%%%%%%%%%%%%%% parameter 1 %%%%%%%%%%%%%%%%%%%%%%%%
SamplingRate = 20000.;

%%%%%%%%%%%%%%%%%%%%% call data %%%%%%%%%%%%%%%%%%%%%%%%%

file = strcat(directory,filename)
filenameold = filename ;
%%% determine the channel length
data = mreadMCD(file, MEA_layout(1));
if ~exist('t1','var')
    t1=0;
end
if ~exist('t2','var')
    t2=0 ;
end
n_orgData = length(data) ;
data = data(t1 * SamplingRate + 1 : n_orgData - t2 * SamplingRate);

%%%%%%%%%%%%%%%%%%%%%% parameters 1 %%%%%%%%%%%%%%%%%%%%%%%
n_rows=length(data);
DataTime = n_rows / SamplingRate ;

%%%%%%%%%%%%%%%%%%%% read mcd data channel by channel  %%%%%%%%%%%%%%%%%
for channel_index = 1 : 60
    tic
    channel_index
    data = mreadMCD(file, MEA_layout(channel_index)); 
    data = data(t1 * SamplingRate + 1 : n_orgData-t2 * SamplingRate);
    dataBaseLine = smooth(data,1000) ;
    data = data - dataBaseLine ;
    %%% peak detecttion
    mean_data = mean(data) ;
    std_data = std(data) ;

    p = peakfinder(data, 0, -1*std_data*5,-1);
%     plot(data) ,hold on, plot(p,std_data*-5,'o'), hold off 
    SpikeTime{channel_index} = p/SamplingRate ;
    toc
end

filename = ['MEA_spike_' filename(1:length(filename)-4) '.mat']
file = strcat(directory,filename) 
save(file,'SpikeTime')

%%%%%%%%%%%%%%%%%%%%%%% analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%
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
imagesc(t_line,[1:60],~h)
colormap('gray')
% F1=getframe(gcf) ; get the image from figure 'getframe(gcf)' or 'getframe'
% imwrite(F1.cdata,file,'jpg')
saveas(gca,file,'jpg'); % get the image from figure 'gca' or 'gcf'
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
colormap('gray')
% F1=getframe ;
% saveas(F1.cdata,file,'jpg')
saveas(gca,file,'jpg');
'finish burst binning'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% histogram analysis,  sum spikes along the channel axis ,seperate burst by BurstDetect_Tina01
name = ['BurstHist_' filename(10:length(filename)-4)]
file = [strcat(directory,name) '.tif' ] 

SumBurst = sum(BinningSpike) ;
SumBurstCutInd = find(SumBurst < 2.5) ;
SumBurstCut = SumBurst ;
SumBurstCut(SumBurstCutInd) = 0 ;
SumBurstTime = t_line(SumBurstCut > 0.5) ;
set_SpikeInterval = 0.2 ;
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
    if abs(b(1) - b(length(b))) > 0.8
        SepBurst(i).his = SumBurst(bb-20:bbb+20) ;
        SepBurst(i).time = t_line(bb-20:bbb+20) ;
    else
        SepBurst(i).his = SumBurst(bb-20:bb+180) ;
        SepBurst(i).time = t_line(bb-20:bb+180) ;
    end
    bar(SepBurst(i).time,SepBurst(i).his) ;
    set(gcf,'color',[1 1 1]) % set gackground as white
    F1=getframe(gcf) ; 
    imwrite(F1.cdata,file,'writemode','append')
end
name = ['Histo_' filename(10:length(filename)-4) '.mat']
if ~exist('Histo', 'dir')
    mkdir(directory,'Histo') 
end;
DirectorySub = [directory 'Histo\'] ;
file = [strcat(DirectorySub,name) '.tif' ] 
save(file,'SepBurst')

