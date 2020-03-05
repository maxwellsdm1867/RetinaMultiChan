HistNum = [] ;
for i = 1 : 23
    temp = StiMat(:,i) ;
    HistNum = [HistNum,temp'] ;
end

p = peakfinder(HistNum,5,10,1) ;
[burst , BurstTime] = BurstDetect_Tina01(p,200,100) ;

n = length(burst) ;
for i = 1 : n
    a = burst{i} ;
    s(i) = length(a) ;
    ind0 = find (p == a(1)) ;
    ind = find (p == a(s(i))) ;
    b0(i) = p(ind0) ;
    b1(i) = p(ind+1) ;
    b2(i) = p(ind+2) ;
    b3(i) = p(ind+3) ;
    b4(i) = p(ind+4) ;
end

s = 6000 ;
StiMat3Bur = zeros(20, s+1) ;
t = [1:s+1]/200 ;

for i = 1 : 20
%     temp = zeros(1,max(b4-b0)+1) ;
%     temp2 = Hist(b2(i):b4(i)) ;
%     temp(1:length(temp2)) = temp2 ;
    StiMat3Bur(i,:) = Hist(b1(i):b1(i)+s) ;
    subplot(4,5,i)
    plot(t,StiMat3Bur(i,:))
    title(num2str(i))
end
