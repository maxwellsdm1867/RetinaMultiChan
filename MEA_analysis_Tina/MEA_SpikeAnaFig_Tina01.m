close all
clear all
clc

cd('E:\Tina\data_MEA\BMI\G9-1\Histo\BurstAna') ; % the folder of the files
all_file = dir('*.mat') ; % change the type of the files which you want to select.
n_file = length(all_file) ; 

% cd('E:\Tina\Dropbox\MatlabCode\biocam')
for k = 1 : n_file
%     clearvars -except all_file n_file k
    file = all_file(k).name ;
    load(file)
%     [pathstr, name, ext] = fileparts(file)
%     directory = [pathstr,'\']
%     filename = [name,ext]
%     
    
   n_burst = length(BurstAna) ;
    for i = 2 : n_burst
        NB(i) = BurstAna(i).NB ;
        TauB(i) = BurstAna(i).TauB ;
        TauBr(i) = BurstAna(i).TauBr ;
    end
    MeanNB(k) = mean(NB) ;
    StdNB(k) = std(NB) ;
    MeanTauB(k) = mean(TauB) ;
    StdTauB(k) = std(TauB) ; 
end
figure
errorbar(MeanNB,StdNB)
figure
errorbar(MeanTauB,StdTauB)

'finish MEA_SpikeTimeAnalysis_Folderuse';