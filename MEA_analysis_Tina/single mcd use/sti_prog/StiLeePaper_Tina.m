% this program is for Kyoung J. Lee' paper PRL 108, 138103 (2012)
% by Tina, 2012 Oct 20

close
clear
% I_sti = 0 ; % which time series u dont want? start time
% F_sti = 400 ; % which time series u dont want? finish time

% [filename,directory] = uigetfile('*.mcd','Pick the AXON file'); 
% file = strcat(directory,filename)
file = 'E:\Tina\data_MEA\sti_lee\20121025sti\20121005(MEA(G12-9))\20121005(MEA(G12-9))_100_20121025001.mcd'
data = -mreadMCD(file, MEA_layout(1));
% data(I_sti*SamplingRate+1 : F_sti*SamplingRate) = [] ;
SamplingRate = 20000.;
n_rows=length(data);

t_data = [1:n_rows] / SamplingRate ;
DataTime = n_rows / SamplingRate ;

ProbInt = 30; % probing stimulation frequence =0.2Hz , interval = 5sec.
probIntPoint = SamplingRate * ProbInt ;
data_I = 4 ; % (sec) data after stimulation 
data_F = 1 ; % (sec) data before stimulation

[pks,locs] = findpeaks(data, 'MINPEAKHEIGHT',4E-4,'MINPEAKDISTANCE',probIntPoint/2) ;
%plot(t_data, data , locs/SamplingRate ,2E-3 , 'o') 

n_sti = length(locs) ; % how many stimulations?
IStiI =  floor(min(locs(2:n_sti)-locs(1:n_sti-1))) ; % inter sti interval
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
        startP = locs(i)+data_I*SamplingRate ;
        a = data(startP : startP+SepDataL) ;
        b = smooth(data(startP : startP+SepDataL),1000);
        c = (a-b)' ;
        ReBurst_one(i,:) = a ;
        p = peakfinder(c, 0, std(c)*-4,-1);
        [n,xout] = hist(p,BinningTime);
        Bin_ReBurst_one(i,:) = n ;
    end
    toc
end
ReBurst_all = ReBurst_all + ReBurst_one ;
Bin_ReBurst_all = Bin_ReBurst_all + Bin_ReBurst_one ;

ss = size(ReBurst_all) ;
ss = ss(2);
tt = [1:ss]/20000 ;
figure
imagesc(tt,[1:n_sti],ReBurst_all)

ss = size(Bin_ReBurst_all) ;
ss = ss(2);
tt = [1:ss]/200 ;
figure
imagesc(tt,[1:n_sti],Bin_ReBurst_all)
'end StiLeePaper_Tina'
