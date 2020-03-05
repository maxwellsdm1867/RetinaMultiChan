%%% SubburstFiringLoc
%%% run after biocam_SpikeTimeAnalysis_Tina02.m 
%%% data save as BurstAna(i).TauB, BurstAna(i).TauBr, BurstAna(i).NB,
%%%              BurstAna(i).TauSb, BurstAna(i).TauSbR, BurstAna(i).NSb
%%% Jan 2013, By Tina
% close all
% clear all
%%%%%%%%%%%%%%% parameter %%%%%%%%%%%%%%%%%%
% [filename,directory] = uigetfile('*.mat','Pick a BurstHist_xxx.mat file'); 
file = strcat(directory,filename)
load(file)
% 
% [filename,directory] = uigetfile('*.mat','Pick a BinningSpike_xxx.mat file'); 
% file = strcat(directory,filename)
% load(file)
SepBurst(1) = [] ;
if ~exist('BurstAna', 'dir')
    mkdir(directory,'BurstAna') ;
end
    DirectorySub = [directory 'BurstAna\'] ;

DataTime = 300 ;
t_line = [0:0.005:DataTime] ;
temp = 0 ;

%%%%%%%%%%%%%%%%%%%%%%%%%% burst ana %%%%%%%%%%%%%%%%%%%%%%%%%
n_burst = length(SepBurst);
for k = 1 : n_burst
close all
clearvars -except BinningSpike SepBurst directory filename DirectorySub...
           t_line DataTime k BurstI BurstF BurstAna temp all_file n_file m
   
    Time = SepBurst(k).time;
    Data = SepBurst(k).his ; %%% t the first burst only
    Time(Data<3) = [] ;
    Data(Data<3) = [] ;
    
    b = smooth(Data,10) ;
    c = smooth(b) ;
    [SubburstAmp,SubburstLoc] = findpeaks(c,'MINPEAKHEIGHT',10,'MINPEAKDISTANCE',10) ; % subburst location and subburst amplitude
    plot (Data)
    hold on
    plot(SubburstLoc,SubburstAmp,'or')
    hold off
    
    if isempty(SubburstAmp)
        [SubburstAmp,SubburstLoc] = max(c) ;
    
    end
    BurstI(k) =  Time(1);
    BurstF(k) = Time(length(Time)) ;
    BurstAna(k).NB = sum(Data) ;
    BurstAna(k).SubT = Time(SubburstLoc);      
    BurstAna(k).TauB = Time(length(Time)) - Time(1) ;
    BurstAna(k).TauBr = Time(1) - temp ;
    temp =  Time(length(Time)) ;
    %%%%%%%%%%%%%%% detect subburst %%%%%%%%%%%%%%%%%%%%
    if length(SubburstAmp) < 1.5
        SubburstLoc = [1 length(c)];
    else
        for i = 1 : length(SubburstAmp)-1
            [minC,minCL] =min(c(SubburstLoc(i):SubburstLoc(i+1))) ;
            minCLoc(i) = minCL +SubburstLoc(i) ;
        end       
        SubburstLoc = [1 minCLoc length(c)];
    end
% 
    n_burstburst = length(SubburstAmp) ;

    %%%%%%%%%%%%%%%%%% subburst ana %%%%%%%%%%%%%%%%%%

    for i = 1 : n_burstburst

        Time1 = SubburstLoc(i) ;
        Time2 = SubburstLoc(i+1) ;
        
%         t1Ind = find(t_line == Time(Time1)) ;
%         t2Ind = find(t_line == Time(Time2)) ;
%         BurstMat = BinningSpike(:,t1Ind:t2Ind) ;
%         d = smooth(sum(BurstMat)) ;
%         
%         dt = 1:length(d) ;
% %         f = fit(dt',d,'gauss1') ; %% fit gaussian curve
%         try 
%         f = fit(dt',d,'gauss1') ;
%         catch
%         f = fit(dt',d,'gauss2') ;
%         end
%           
%         Time3 = fix(f.b1-2*f.c1) ;
%         if Time3 < 1 ;
%             Time3 = 1;
%         end
%         Time4 = fix(f.b1+2*f.c1) ;
%         if Time4 > length(d) ;
%            Time4 = length(d) ;
%         end
% 
%         BurstMat = BurstMat(:,Time3:Time4) ;
        dump3(i) = sum(Data(Time1:Time2)) ;
%         dump4(i) = 2 * f.c1 ;
    end
BurstAna(k).NSb = dump3 ;
% BurstAna(k).SubWidth = dump4 ;
BurstAna(k).SubWidth = diff(SubburstLoc) ;
end
BurstAna(1) = [] ;

name = ['BurstAna_' filename(6:length(filename)-4) '.mat'];
file = [DirectorySub name]
save(file,'BurstAna')

'end biocam_SpikeAna_Tina01'