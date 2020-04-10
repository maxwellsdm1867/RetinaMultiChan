close all 
clear all
clc

merge_file_dir = 'D:\sfn2019\contrast\190612_contrast_merge.mat';
rslt_dir = 'D:\sfn2019\contrast\190612';
mcd_dir = 'D:\data\190612\';
name_seq = {'tr1c02n0','tr1c03n0','tr1c04n0','tr1c05n0','tr1c01n0','tr2c02n4','tr2c02n4_2','sta'};
trial_time = [215.4 215.4 215.4 215.4 215.4 215.4 215.4 305.4];
sort_extract(merge_file_dir,name_seq,trial_time,rslt_dir,mcd_dir)