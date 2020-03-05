%%% laod files 

close all
clear all
clc

cd('E:\MEA result\Webber_2014\20140411(MEA(F5-18))') ; % the folder of the files
all_file = subdir('*.mcd') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file) ; 

cd('E:\Tina\Dropbox\MatlabCode\MEA\MEA_analysis_Tina')
for m = 1 : n_file
    clearvars -except all_file n_file m
    file = all_file(m).name ;
    [pathstr, name, ext] = fileparts(file)
    directory = [pathstr,'\']
    filename = [name,ext]
    MEA_SpikeTimeAnalysis_Tina02  % the name of program
end
'finish MEA_SpikeTimeAnalysis_Folderuse';