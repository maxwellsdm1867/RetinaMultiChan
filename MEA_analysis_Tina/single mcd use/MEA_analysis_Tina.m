% analysis MEA data wehn we have m file, which is create by 'MEA_SpikeTimaRecord_Tina'
% by Tina 2012, May

clc
clear all
close all
[filename,directory] = uigetfile('*.mat','Pick a MEA_spikeXXX.. file'); 
file = strcat(directory,filename) ;
load(file) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% binning data with 5ms interval, use function 'BinningData_Tina01'
DataTime = input('how many seconds of this data?  ')
%DataTime = 600 ;
t_line = [0:0.005:DataTime] ;
[BinningSpike] = BinningData_Tina01(SpikeTime,DataTime) ;
imagesc(BinningSpike)

name = ['BinningSpike_' filename(10:length(filename)-4)]
file = strcat(directory,name) 
save(file,'BinningSpike')
file = [file '.jpg'] 
h = mat2gray(BinningSpike) ;
imagesc([0:0.005:DataTime],[1:60],~h)
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
    [burst , BurstTime] = BurstDetect_Tina01(aa,set_SpikeInterval,set_SpikeInBurst) ;
    MEABurst{i} = burst ;
    MEABurstTime{i} = BurstTime ;
end
name = ['MEA_Burst_' filename(10:length(filename)-4)]
file = strcat(directory,name)
save(file,'MEABurstTime')
'finish burst detect'

[BinningBurst] = BinningData_Tina01(MEABurstTime,DataTime) ;
imagesc(BinningBurst)
name = ['BinningBurst_' filename(10:length(filename)-4)]
file = strcat(directory,name)
save(file,'BinningBurst')

file = [file '.jpg'] ;
h = mat2gray(BinningBurst) ;
imagesc([0:0.005:t_line],[1:60],~h)
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
for i = 2 :length(HistBurst)-1
    b = HistBurst{i} ;
    bb = find(t_line == b(1)) ;
    bl = length(b) ;
    bbb = find (t_line == b(bl)) ;
    if abs(b(1) - b(length(b))) > 1.0
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
name = ['Histo_' filename(10:length(filename)-4)]
file = strcat(directory,name)
save(file,'SepBurst')

'finish burst histgram'
