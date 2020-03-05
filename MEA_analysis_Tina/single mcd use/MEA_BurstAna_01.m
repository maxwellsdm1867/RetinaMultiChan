%%% MEA IBI histogram
clear all
close all

[filename,directory] = uigetfile('*.mat','Pick a Histo_xxx.mat file'); 
file = strcat(directory,filename)
load(file)

n_burst = length(SepBurst) ;
BTF = 0 ;
BTI = 0 ;

for i = 2 : n_burst

%%%%%%%%%%%%%%%%%%%%%% IBI %%%%%%%%%%%%%%%%%%%%%%%%%
    Data = SepBurst(i).his ;
    Time = SepBurst(i).time ;
    BurstNSpk(i) = sum(Data) ;
          
    BTI = Time(1) ; % srart of a burst
    BurstStart(i) = BTI ;
    IBI(i) = BTI - BTF ;
    
    BTF = Time(length(Time)) ; % the end of a burst
    BurstStop(i) = BTF;
    
end
IBI(1) = [] ;

x =  IBI ; 
y = BurstNSpk(2:n_burst) ; 
DumpInd = find(IBI<1) ;
x(DumpInd) = [] ;
y(DumpInd) = [] ;
% name = ['BurstAna_' filename(6:length(filename)-4) '.mat']
% file = strcat(directory,name)
% save(file,'SepBurstSta')


%%%%%%%%%%%%%%%%%%% save IBI vs. spikes figure %%%%%%%%%%%%%%%%%%%
name = ['IBI_' filename(7:length(filename)-4)]
if ~exist('IBI', 'dir')
    mkdir(directory,'IBI') ;
end
DirectorySub = [directory 'IBI\'] ;
file = [DirectorySub name '.tif'] ;
% plot(x,y,'o')
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

%%%%%%%%%%%%%% save distribution of spikes in a burst %%%%%%%%%%%%%%%%%%
name = ['NSpikes_' filename(7:length(filename)-4)]
DirectorySub = [directory 'IBI\'] ;
file = [DirectorySub name '.tif'] ;
[n,xout] = hist(y);
bar(xout,n)
set(gcf,'color',[1 1 1]) % set gackground as white
F1=getframe(gcf) ; 
imwrite(F1.cdata,file,'writemode','append')

%%%%%%%%%%%%%% save distribution of IBI in a burst %%%%%%%%%%%%%%%%%%
name = ['NIBI_' filename(7:length(filename)-4)]
DirectorySub = [directory 'IBI\'] ;
file = [DirectorySub name '.tif'] ;
[n,xout] = hist(x);
bar(xout,n)
set(gcf,'color',[1 1 1]) % set gackground as white
F1=getframe(gcf) ; 
imwrite(F1.cdata,file,'writemode','append')

'end'