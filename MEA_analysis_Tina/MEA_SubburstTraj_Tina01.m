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

k = 18 ;

name = ['Subburst' num2str(k) 'Trace_' filename(10:length(filename)-4)] ;
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

for i = 1 : length(SubburstAmp)-1
    [minC,minCL] =min(c(SubburstLoc(i):SubburstLoc(i+1))) ;
    minCLoc(i) = minCL +SubburstLoc(i) ;
end
SubburstLoc = [1 minCLoc length(c)];
N_Subburst = length(SubburstAmp) ;


for i = 1 : N_Subburst

        tt1 = SubburstLoc(i) ;
        tt2 = SubburstLoc(i+1) ;
        t1Ind = find(t_line == tt(tt1)) ;
        t2Ind = find(t_line == tt(tt2)) ;
        BurstMat = BinningSpike(:,t1Ind:t2Ind) ;

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
        quiver(CenterX(1:length(CenterX)-1),CenterY(1:length(CenterY)-1),diff(CenterX),diff(CenterY),0.5)
        % plot line with different color over time
        z = zeros(size(colorlength));
        tcolor = colorlength;  % This is the color, vary with x in this case.     
        surface([CenterX;CenterX],[CenterY;CenterY],[z;z],[tcolor;tcolor],'edgecol','interp','linew',2);
        hold off
        axis([1, 8, 1, 8])
        set(gca,'YDir','reverse')

        tNow = num2str(tt(tt1),'%7.3f') ;
        SubburstID = num2str(i) ;
        tit = [filename(13:length(filename)-10) ', SubburstID=' SubburstID ', t=' tNow];
        title(tit)
        set(gcf,'color',[1 1 1]) % set gackground as white
        F1=getframe(gcf) ; 
        imwrite(F1.cdata,file,'writemode','append')

        %%%%% image
        N_fig = tt2-tt1+1 ;
        SpkMat = zeros(8) ;
        scale = fix(255/N_fig) ;
        for j = 1 : N_fig
            SpkNow = BurstMat(:,j);
            SpkInd = find(SpkNow) ;
            SpkMatInd = sub2ind(size(SpkMat),ChY(SpkInd),ChX(SpkInd)) ; % convert the (x,y) to linear index
            SpkMatAmp = j*scale ;
            SpkMat(SpkMatInd) = SpkMatAmp; % mark the firing
            image(SpkMat),colormap(jet(256)),colorbar
            tNow = num2str(tt(tt1+j-1),'%7.3f') ;
            SubburstID = num2str(i) ;
            tit = [filename(13:length(filename)) ', SubburstID=' SubburstID ', t=' tNow];
            title(tit)
            set(gcf,'color',[1 1 1]) % set gackground as white
            F1(j)=getframe(gcf) ; 
            imwrite(F1(j).cdata,file,'writemode','append')
        end
end



'end MEA_subburstTraj'