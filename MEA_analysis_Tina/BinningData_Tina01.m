%%% binning data by Tina 2012, May

function [BinningSpike] = BinningData_Tina01(SpikeTime,DataTime) 

BinningInterval = 5/1000 ; % (s) binning duration = 5ms
BinningTime = [0:BinningInterval:DataTime] ;% # 0~600s with 5ms interval
BinningSpike = zeros(60,length(BinningTime)) ;

for i = 1:60
    [n,xout] = hist(SpikeTime{i},BinningTime) ;
    BinningSpike(i,:) = n ;
end

