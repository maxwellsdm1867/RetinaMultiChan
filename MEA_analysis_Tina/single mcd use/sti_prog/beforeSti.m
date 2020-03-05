%%%time deries before probing and training, 300s => 30s *10
close all
clear all

file = 'E:\Tina\data_MEA\sti_lee\20121029sti\20121012(MEA(F10_10))_100_-30_20121029001.mcd'
data1 = -mreadMCD(file, MEA_layout(1));
SamplingRate = 20000.;
data1(300*20000+1 : length(data1)) = [] ;
tt = [1:length(data1)];

for channel_index = 1 : 30
    channel_index
    tic
    data = mreadMCD(file, MEA_layout(channel_index));
    data(300*20000+1 : length(data)) = [] ;
    data1 = data1 + data ;
end


TraceBefore = zeros(10,SamplingRate*30);

for i = 1 :10
    TraceBefore(i,:) = data1((i-1)*30*SamplingRate+1 : i*SamplingRate*30);
end
% TraceBefore2 = TraceBefore ;
% overser = TraceBefore > 1E-4;
% TraceBefore(overser) = 1E-4 ;
% overser = TraceBefore < -1E-4;
% TraceBefore(overser) = -1E-4 ;
figure
imagesc(tt ,[1:10],TraceBefore)