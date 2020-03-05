close all 
clear all
clc

merge_file_dir = 'D:\analyse\190529\varability_merge.mat';
rslt_dir = 'D:\analyse\190529\var_test';
mcd_dir = 'D:\data\190529\';
name_seq = {'onoff_start','tr2c02n3','tr2c02n3-2','tr2c03n3','tr2c03n3-2','sta'};
trial_time = [90.2 215.4 215.4 215.4 215.4 305.4];
sort_extract(merge_file_dir,name_seq,trial_time,rslt_dir,mcd_dir)