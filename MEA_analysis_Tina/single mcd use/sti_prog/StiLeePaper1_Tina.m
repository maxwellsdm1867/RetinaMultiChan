% this program is for Kyoung J. Lee' paper PRL 108, 138103 (2012)
% by Tina, 2012 Oct 20

close
clear
% I_sti = 0 ; % which time series u dont want? start time
% F_sti = 400 ; % which time series u dont want? finish time

% [filename,directory] = uigetfile('*.mcd','Pick the AXON file'); 
% file = strcat(directory,filename)
file = 'E:\Tina\data_MEA\sti_lee\20121031sti\20121005(MEA(A19-18))_Tina_20121031001\20121005(MEA(A19-18))_Tina_20121031004.mcd'
data = -mreadMCD(file, MEA_layout(1));
% data(I_sti*SamplingRate+1 : F_sti*SamplingRate) = [] ;
SamplingRate = 20000.;
n_rows=length(data);

t_data = [1:n_rows] / SamplingRate ;
DataTime = n_rows / SamplingRate ;

ProbInt = 30; % probing stimulation frequence =0.2Hz , interval = 5sec.
probIntPoint = SamplingRate * ProbInt ;
data_I = 0 ; % (sec) data after stimulation 
data_F = 0 ; % (sec) data before stimulation

[pks,locs] = findpeaks(data, 'MINPEAKHEIGHT',4E-4,'MINPEAKDISTANCE',probIntPoint/2) ;
%plot(t_data, data , locs/SamplingRate ,2E-3 , 'o') 

n_sti = length(locs) ; % how many stimulations?
IStiI =  min(locs(2:n_sti)-locs(1:n_sti-1)) ; % inter sti interval
SepDataL = IStiI - (data_I+data_F)*SamplingRate ; % one sti data length , ignore the first 3 sec after sti and 1 sec before sti
ReBurst_one = zeros(n_sti,SepDataL+1);
ReBurst_all = ReBurst_one ;

BinningTime = [0:100:SepDataL] ; % binning with datapoint (20K Hz to 200Hz (5ms))

Bin_ReBurst_one = zeros(n_sti,length(BinningTime)) ; %binning spikes
Bin_ReBurst_all = Bin_ReBurst_one ;


for channel_index = 1 : 60
    channel_index
    tic
    data = mreadMCD(file, MEA_layout(channel_index));
%     data(I_sti*SamplingRate+1 : F_sti*SamplingRate) = [] ;
    
    for i = 1 : n_sti
        startP = locs(i)+data_I*SamplingRate+1 ;
        a = data(startP : startP+SepDataL) ;
        b = smooth(a,1000);
        c = (a-b)' ;
        ReBurst_one(i,:) = a ;
        c(1:1*50000)=[] ;
        p = peakfinder(c, 0, std(c)*-4,-1)+50000 ;
        [n,xout] = hist(p,BinningTime);
        Bin_ReBurst_one(i,:) = n ;
    end
    toc
    ReBurst_all = ReBurst_all + ReBurst_one ;
    Bin_ReBurst_all = Bin_ReBurst_all + Bin_ReBurst_one ;
end


ss = size(ReBurst_all) ;
ss = ss(2);
tt = [1:ss]/20000 ;
figure
imagesc(tt,[1:n_sti],ReBurst_all)
SclReBurst_all = ReBurst_all ;
overser = ReBurst_all > 1E-4;
SclReBurst_all(overser) = 1E-4 ;
overser = ReBurst_all < -1E-4;
SclReBurst_all(overser) = -1E-4 ;
figure
imagesc(tt,[1:n_sti],SclReBurst_all)
% 
% name = ['MEA_spike_' filename(1:length(filename)-4) '.mat']
% file = strcat(directory,name) 
% save(file,'ReBurst_all')

ss = size(Bin_ReBurst_all) ;
ss = ss(2);
tt = [1:ss]/200 ;
figure
imagesc(tt,[1:n_sti],Bin_ReBurst_all)

% name = ['MEA_spike_' filename(1:length(filename)-4) '.mat']
% file = strcat(directory,name) 
% save(file,'Bin_ReBurst_all')
'end StiLeePaper_Tina'
