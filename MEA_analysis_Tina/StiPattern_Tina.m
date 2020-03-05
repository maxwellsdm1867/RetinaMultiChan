%%% plot the burst time for a stimulation experiment
%%% read BinningSpike__XXXXXX.mat
%%% by Tina , Jun 2012

clear all
close all
clc

BeforeSti = 400 ; % (data point / rate ) un-stimulation time
StiInt = 240 ; % stimulation interval 
Rate = 200 ; % sampling rate

cd('E:\Tina\data_MEA\sti_20120619\BinningSpike') ;
all_file = dir('*.mat') ;
n_file = length(all_file) ; 
AppSpike = [] ;
sss = 0 ;

for i = 1 : n_file
    file = all_file(i).name ;
    [pathstr, name, ext] = fileparts(file) ;
    directory = [pathstr,'\'] ;
    filename = [name,ext] 
    load(filename) ;
    a{i} = BinningSpike ;
    SumSpike = sum(BinningSpike) ;
    ss(i) = length(SumSpike) ;
    sss = sss+length(SumSpike) ;
    AppSpike = [AppSpike , SumSpike] ;
end
nx = StiInt*Rate ; %(time * rate)
ny = ceil(length(AppSpike) / nx) ;
ns = nx*ny ;
AppSpike3 = AppSpike((BeforeSti-StiInt)* Rate : length(AppSpike)) ;
dump = zeros(1,ns - length(AppSpike3)) ;
AppSpike2 = [AppSpike3 , dump ] ;
Th = 10 ;
IndTh = find(AppSpike2 < 10) ;
AppSpike2(IndTh) = 0 ;
% IndTh2 = find(AppSpike2 > 10) ;
% AppSpike2(IndTh) = 250 ;
StiMat = reshape(AppSpike2 , nx ,ny) ;
StiMatInd = StiMat > 10 ;
h = mat2gray(StiMatInd') ;
t = [0.005:0.005:StiInt] ;
% scrsz = get(0,'ScreenSize');
% hfig=(figure);
% set(hfig,'position',scrsz);
imagesc(t,[1:ny],~h)
colormap('gray')
xlabel('time/s')
ylabel('trails')

directory = cd ;
name = ['StiPattern',filename(13:length(filename)-8)];
file = [directory ,'\',name '.jpg' ] 
saveas(gca,file,'jpg') ;
save(name,'StiMat');

for i = 1:23
subplot(5,5,i)
plot(t,StiMat(:,i))
title(num2str(i))
% xlabel('time/s')
% ylabel('# of spikes')
end

name = ['StiHistogram',filename(13:length(filename)-8)];
file = [directory ,'\',name '.jpg' ] 
saveas(gca,name,'jpg') ;


%%%%%%%%%%% paper result 

p = peakfinder(AppSpike2,5,10,1) ;
[burst , BurstTime] = BurstDetect_Tina01(p,200,100) ;

n = length(burst) ;
for i = 1 : n
    aa = burst{i} ;
    s(i) = length(aa) ;
    ind0 = find (p == aa(1)) ;
    ind = find (p == aa(s(i))) ;
    b0(i) = p(ind0) ;
    b1(i) = p(ind+1) ;
    b2(i) = p(ind+2) ;
    b3(i) = p(ind+3) ;
    b4(i) = p(ind+4) ;
end

s = max(b4-b1) ;
StiMat3Bur = zeros(20, s+1) ;
t = [1:s+1]/200 ;

for i = 1 : 20
%     temp = zeros(1,max(b4-b0)+1) ;
%     temp2 = Hist(b2(i):b4(i)) ;
%     temp(1:length(temp2)) = temp2 ;
    StiMat3Bur(i,:) = AppSpike2(b1(i)+2:b1(i)+s) ;
    subplot(4,5,i)
    plot(t,StiMat3Bur(i,:))
    title(num2str(i))
end

name = ['StiHistogram2',filename(13:length(filename)-8)];
file = [directory ,'\',name '.jpg' ] 
saveas(gca,name,'jpg') ;

StiMatInd = StiMat3Bur > 5 ;
h = mat2gray(StiMat3Bur) ;
t = [0.005:0.005:length(StiMat3Bur(1,:))] ;
% scrsz = get(0,'ScreenSize');
% hfig=(figure);
% set(hfig,'position',scrsz);
imagesc(t,[1:ny],~h)
colormap('gray')
xlabel('time/s')
ylabel('trails')

'finish StiPattern'