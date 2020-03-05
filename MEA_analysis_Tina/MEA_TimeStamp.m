%%% timestamp from MC_rack
close all
clear all
clc

cd('E:\Tina\data_MEA\BMI\G9-1\timestamp') ; % the folder of the files
all_file = subdir('*.dat') ; % change the type of the files which you want to select.
n_file = length(all_file) ; 

for i = 1 : n_file
    file = all_file(i).name ;
    A = importdata(file);
    SpikeTime{i} = A.data ;
end

file = [cd '\SpikeTime_100.12.16(MEA(G9-1))_20120109004.mat'];
save(file,'SpikeTime')

'finish SpikeTime';

