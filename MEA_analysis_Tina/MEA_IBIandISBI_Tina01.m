%%% MEA IBI(inter burst interval) and  ISBI(inter subburst interval) distribution
clear all
close all

[filename,directory] = uigetfile('*.mat','Pick a Histo_xxx.mat file'); 
file = strcat(directory,filename)
load(file)

n_burst = length(SepBurst) ;
BTF = 0 ;
BTI = 0 ;
figure
for i = 2 : n_burst
    %%%%%%%%%%%%%%%%%%%% ISBI %%%%%%%%%%%%%%%%%%%%%
    clearvars a b d f dump
    Data = SepBurst(i).his ;
    Time = SepBurst(i).time ;
    a = smooth(Data,10) ;
    b = smooth(a) ;
    
    plot(SepBurst(i).time,b)
    [pks,locs] = findpeaks(b,'MINPEAKHEIGHT',10,'MINPEAKDISTANCE',10);
       
    TimePks = Time(locs) ;
    DataPks = Data(locs) ;
    SepBurstSta(i).time = TimePks ; % subburst peak time
    SepBurstSta(i).pks = DataPks ;  % subburst peak amplitude
    
    n_locs = length(locs) ;
    dump = fix((locs(2:n_locs) + locs(1:n_locs-1))/2) ;
    d = [1, dump' , length(Data)];
    SepBurstSta(i).regi = d'; % subburst region
        
    ISBI = diff(TimePks);
    SepBurstSta(i).ISBI = ISBI ; % subburst inter burst interval
    
    for j = 1 : n_locs
         f(j) = sum(Data(d(j):d(j+1))) ;
    end
    SepBurstSta(i).NSpk = f ; % # of spikes of subburst 
        
    hold on
    plot(TimePks,10,'*')
    hold off
%     pause(0.5)  


%%%%%%%%%%%%%%%%%%%%%% IBI %%%%%%%%%%%%%%%%%%%%%%%%%
    BurstNSpk(i) = sum(Data) ;
          
    BTI = Time(1) ; % srart of a burst
    BurstStart(i) = BTI ;
    IBI(i) = BTI - BTF ;
    
    BTF = Time(length(Time)) ; % the end of a burst
    BurstStop(i) = BTF;
    
end
IBI(1) = [] ;
name = ['IBIAndISBI_' filename(6:length(filename)-4) '.mat']
file = strcat(directory,name)
save(file,'SepBurstSta')

%%%%%%%%%%% save ISBI vs. spikes figure %%%%%%%%%%%%%%%%%%%
close all
name = ['ISBI_' filename(6:length(filename)-4)]
if ~exist('ISBI', 'dir')
    mkdir(directory,'ISBI') 
end;
DirectorySub = [directory 'ISBI\'] ;
file = [strcat(DirectorySub,name) '.tif' ] 
for i = 2 : n_burst
    clearvars e n_e
    e = SepBurstSta(i).NSpk ;
    n_e = length(e) ;
    plot(SepBurstSta(i).ISBI,e(2:n_e),'o')
    tit = ['burst ' num2str(i)];
    title(tit)
    xlabel('ISBI / sec')
    ylabel('# of spikes')
    set(gcf,'color',[1 1 1]) % set gackground as white
    F1=getframe(gcf) ; 
    imwrite(F1.cdata,file,'writemode','append')
%      pause(0.5)
end

%%%%%%%%%%%%%%%%%%% save IBI vs. spikes figure %%%%%%%%%%%%%%%%%%%
name = ['IBI_' filename(6:length(filename)-4)]
if ~exist('IBI', 'dir')
    mkdir(directory,'IBI') ;
end
DirectorySub = [directory 'IBI\'] ;
file = [DirectorySub name '.tif'] ;
x =  IBI ; 
y = BurstNSpk(2:n_burst) ; 
DumpInd = find(IBI<1) ;
x(DumpInd) = [] ;
y(DumpInd) = [] ;
plot(x,y,'o')
[p,S] = polyfit(x,y,1);
[linearfit,delta] = polyval(p,x,S);
FitError = mean(delta) ;
plot(x,y,'s', x,linearfit,'r-')
title(name)
xlabel('IBI / sec')
ylabel('# of spikes')
func = ['y=' num2str(fix(p(1))) 'x+' num2str(fix(p(2))) ', MeanSD=' num2str(fix(FitError))] ; 
text(min(x),max(y),func)
set(gcf,'color',[1 1 1]) % set gackground as white
F1=getframe(gcf) ; 
imwrite(F1.cdata,file,'writemode','append')

'end'