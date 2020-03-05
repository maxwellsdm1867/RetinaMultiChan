%%% plot the spikes in some range between 'SelMin' and 'SelMax'
%%% input a BinningSpike matrix
%%% 
%%% [SelData,t] = MEA_VewRangeSpike(BinningSpike,56,58) ;
%%% imagesc(t,[1:60],~SelData)
%%% set(gcf,'color',[1 1 1])
%%% colormap('gray')


% clear all
% close all
% 
% [filename,directory] = uigetfile('*.mat','Pick a BinningSpike_.. file'); 
% file = strcat(directory,filename) ;
% load(file) ;

function [SelData,t] = MEA_VewRangeSpike(BinningSpike,SelMin,SelMax)

% SelMin = 56 ; % select interval (s)
% SelMax = 58 ; % select interval (s)
t_line = [0:0.005:600];
SelMaxInd = find(t_line ==SelMax) ;
SelMinInd = find(t_line ==SelMin) ;
t = t_line(SelMinInd:SelMaxInd) ;

SelData = BinningSpike(:,SelMinInd:SelMaxInd) ;
% imagesc(t,[1:60],SelData)
% set(gcf,'color',[1 1 1])
% colormap('gray')
%