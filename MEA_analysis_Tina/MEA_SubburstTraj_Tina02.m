%%% MEA subburst trajectory
%%% 27 May 2013, Tina
clear all
close all


[filename,directory] = uigetfile('*.mat','Pick a Histo_xxx.mat file'); 
file = strcat(directory,filename)
load(file)

[filename,directory] = uigetfile('*.mat','Pick a BinningSpike_xxx.mat file'); 
file = strcat(directory,filename)
load(file)

ChX = [ repmat(1,1,6) , repmat(2,1,8) ,repmat(3,1,8) ,repmat(4,1,8) ,repmat(5,1,8) ,repmat(6,1,8) ,repmat(7,1,8) ,repmat(8,1,6)];
ChY = [[2:7], [1:8], [1:8], [1:8], [1:8], [1:8], [1:8], [2:7]];
DataTime = 600 ;
t_line = [0:0.005:DataTime] ;



%%%%%%%%%%%%%% subburst trace %%%%%%%%%%%%%%%
%%%%% create folder
n_sub = length(SepBurst);
if n_sub >10
    n_sub = 10 ;
end

for k = 2 : n_sub
close all
clearvars -except BinningSpike ChX ChY SepBurst directory file filename t_line k

name = ['Subburst' num2str(k) 'Trace_' filename(14:length(filename)-4)] ;
if ~exist('SubburstTrac', 'dir')
    mkdir(directory,'SubburstTrac') ;
end
DirectorySub = [directory 'SubburstTrac\'] ;
file = [DirectorySub name '.tif'] 

tt = SepBurst(k).time;
tt(1:18)=[] ;
a = SepBurst(k).his ; %%% t the first burst only
a(1:18)=[] ;
c = smooth(a,10) ;
plot(c)
[SubburstAmp,SubburstLoc] = findpeaks(c,'MINPEAKHEIGHT',10,'MINPEAKDISTANCE',10) ;
hold on 
plot(SubburstLoc,SubburstAmp,'*')
hold off
xlabel('time / 5ms')
ylabel('# of spikes')
BurstID = num2str(k) ;
tit = [filename(14:length(filename)-4) ', BurstID=' BurstID];
title(tit)
set(gcf,'color',[1 1 1]) % set gackground as white
F1=getframe(gcf) ; 
imwrite(F1.cdata,file,'writemode','append')

if length(SubburstAmp) < 1.5
        SubburstLoc = [1 length(c)];
    else
        for i = 1 : length(SubburstAmp)-1
            [minC,minCL] =min(c(SubburstLoc(i):SubburstLoc(i+1))) ;
            minCLoc(i) = minCL +SubburstLoc(i) ;
        end       
        SubburstLoc = [1 minCLoc length(c)];
    end

N_Subburst = length(SubburstAmp) ;


for i = 1 : N_Subburst

        tt1 = SubburstLoc(i) ;
        tt2 = SubburstLoc(i+1) ;
%         t1Ind = find(t_line == tt(tt1)) ;
%         t2Ind = find(t_line == tt(tt2)) ;
        t1Ind = tt(tt1)/0.005+1 ;
        t2Ind = tt(tt2)/0.005+1 ;
        BurstMat = BinningSpike(:,t1Ind:t2Ind) ;
        d = smooth(sum(BurstMat)) ;
        
        dt = 1:length(d) ;
%         f = fit(dt',d,'gauss1') ; %% fit gaussian curve
        try 
        f = fit(dt',d,'gauss1') ;
        catch
        f = fit(dt',d,'gauss2') ;
        end
        
        tt3 = fix(f.b1-2*f.c1) ;
        if tt3 < 1
            tt3 = 1;
        end
        tt4 = fix(f.b1+2*f.c1) ;
        if tt4 > length(d)
           tt4 = length(d) ;
        end

        BurstMat = BurstMat(:,tt3:tt4) ;
        dump = sum(BurstMat)<2 ;
        BurstMat(:,dump)=[] ;

        %%%%% trace
        BurstMat = logical(BurstMat) ;
        SpkChX = bsxfun(@times,ChX',BurstMat);
        SpkChY = bsxfun(@times,ChY',BurstMat);
        CenterX = ceil(sum(SpkChX) ./ sum(BurstMat)) ;
        CenterX(isnan(CenterX))=[] ;
        CenterY = ceil(sum(SpkChY) ./ sum(BurstMat)) ;
        CenterY(isnan(CenterY))=[] ; 
        % plot dot with different color over time 
        colorlength = [1 : length(CenterX)] ;
        scatter(CenterX, CenterY,5,colorlength);
        hold on
        % plot arrow
%         quiver(CenterX(1:length(CenterX)-1),CenterY(1:length(CenterY)-1),diff(CenterX),diff(CenterY))
        % plot line with different color over time
        z = zeros(size(colorlength));
        tcolor = colorlength;  % This is the color, vary with x in this case.     
        surface([CenterX;CenterX],[CenterY;CenterY],[z;z],[tcolor;tcolor],'edgecol','interp','linew',2);
        colormap(jet)
        hold off
        axis([1, 8, 1, 8])
        set(gca,'YDir','reverse')

        tNow = num2str(tt(tt1),'%7.3f') ;
        SubburstID = num2str(i) ;
        tit = [filename(14:length(filename)-4) ', SubburstID=' SubburstID ', t=' tNow];
        title(tit)
        set(gcf,'color',[1 1 1]) % set gackground as white
        F1=getframe(gcf) ; 
        imwrite(F1.cdata,file,'writemode','append')

        %%%%% image
        N_fig = 20;
        if length(sum(BurstMat)) < N_fig  
            N_fig = length(sum(BurstMat)) ;
        end
        SpkMat = zeros(8) ;
        for j = 1 : N_fig
            SpkNow = BurstMat(:,j);
            SpkInd = find(SpkNow) ;
            SpkMatInd = sub2ind(size(SpkMat),ChY(SpkInd),ChX(SpkInd)) ; % convert the (x,y) to linear index
            SpkMatAmp = j ;
            SpkMat(SpkMatInd) = SpkMatAmp; % mark the firing
            image(SpkMat,'CDataMapping','scaled'),colormap(hot),colorbar
            tNow = num2str(tt(tt1+j-1),'%7.3f') ;
            SubburstID = num2str(i) ;
            tit = [filename(13:length(filename)) ', SubburstID=' SubburstID ', t=' tNow];
            title(tit)
            set(gcf,'color',[1 1 1]) % set gackground as white
            F1(j)=getframe(gcf) ; 
            imwrite(F1(j).cdata,file,'writemode','append')
        end
end
end


'end MEA_subburstTraj'