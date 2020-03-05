% this is the main file for analyisis
%% load files

close all
clear all
clc
datez = '190612';
%reset all the path, please check
spike_path = ['D:\sfn2019\cut\' datez];
analog_path = ['D:\data\' datez];
result_path = ['D:\sfn2019\analyse\' datez];
seq_path = 'D:\data\190318\seq.mat';
sta_path = ['D:\sfn2019\cut\' datez '\STA.mat'];
STA_ana=['D:\data\' datez '\sta.mat'];
load(seq_path)
analyse.seq = seq; %save seq for later

%conditions ----change every time
contrast =2:3;
noise =0:5;
conditions = length(noise)*length(contrast);
tr= 'tr3';
comp = cell(3,conditions);
s = 1;
%load file name, spike& channel ID, analog data
for i = contrast
    for j = noise
        filename1 = [tr 'c0' num2str(i) 'n' num2str(j)  '.mat'];%spike file
        load([spike_path '\' filename1])
        comp{1,s} = Spikes;
        comp{2,s} = filename1(1:end-4);
        %filename2 = [ 'c0' num2str(i) 'n' num2str(j)  '.mat'];%analog file
        filename2 = [tr 'c0' num2str(i) 'n' num2str(j)  '.mat'];%analog file
        load([analog_path '\' filename2] )
        comp{3,s} = a_data(2,:);
        s = s+1;
    end
end

%% compute the spike trigger average here
[analyse.sta analyse.staC]= arthur_sta(sta_path,300,STA_ana);


%% process one condition at once
for cdt = 1:conditions
    %cdt = 1;%comment this one later
    tmp_spk = comp{1,cdt};
    tmp_name = comp{2,cdt};
    tmp_analog = comp{3,cdt};
    
    analyse.filename = tmp_name;
    analyse.analog = tmp_analog;
    analyse.sampling_rate = 20000;
    analyse.num_unit = size(tmp_spk,2);
    analyse.data_time = length(tmp_analog)/analyse.sampling_rate;% in seconds
    num_unit= size(tmp_spk,2);
    
    
    %id processing
    %first 2 digits is the channel id and the last digit is the unit id
    id_cell = cell(1,num_unit);
    id_array = zeros(1,num_unit);
    for f= 1:num_unit
        id_cell{1,f} = [num2str(tmp_spk{3,f}) num2str(tmp_spk{2,f})];
        id_array(1,f)=str2num( id_cell{1,f});
    end
    analyse.id = id_array;
    
    %spike count and spike time
    sc = zeros(1,num_unit);
    for k = 1:num_unit
        sc(1,k) = length(tmp_spk{1,k});
    end
    analyse.spike_count = sc;
    analyse.spike_time = tmp_spk(1,:);
    
    %binning spike(raster plot)
    DataTime = analyse.data_time;
    BinningInterval = 0.040;%sec as in 25 Hz
    BinningTime = [0 :  BinningInterval : DataTime];
    all_binspk =zeros(num_unit,length(BinningTime) );
    for m = 1:num_unit
        spike_time = tmp_spk{1,m};
        [BinningSpike] = BinSpk1(BinningInterval,spike_time,DataTime);
        all_binspk(m,:) = [BinningSpike];
    end
    analyse.BinningSpike = all_binspk;
    analyse.BinningTime = BinningTime;
    %inter spike interval
    intervals = cell(1, num_unit);
    for n = 1:num_unit
        k1 = diff(tmp_spk{1,n});
        intervals{1,n} = k1;
    end
    analyse.ISI = intervals;
    
    %Time shifted mutual information, peak height, peak location in time
    a_data = analyse.analog;
    sorted_Spikes = analyse.spike_time;
    Infos = analyse.sampling_rate;
    
    [MI_for_all_chan,t] =TSMI_sorted(a_data,sorted_Spikes,Infos,seq,11);
    analyse.TSMI = MI_for_all_chan;
    analyse.mi_t = t;
    peak_heights = zeros(1,num_unit);
    peak_locs = zeros(1,num_unit);%index
    MI_single = zeros(1,num_unit);
    for p = 1:num_unit
        tmp = MI_for_all_chan{1,p};
        [pks,locs] = max(tmp(130:170));
        if ~isempty([pks,locs])
            peak_heights(1,p) = pks;
            peak_locs(1,p) = locs+130;
        end
        MI_single(p) = tmp(find(t==0));
    end
    analyse.miph = peak_heights;
    analyse.mipl = peak_locs;
    analyse.misv = MI_single;
    thrs = 150;
    rpt = 2;
    [pre_step] = cutspk(a_data,thrs,sorted_Spikes,rpt,Infos);
    analyse.pre = pre_step;
    result_mi = spike_entropy(a_data,sorted_Spikes,Infos,seq );
    analyse.spike_entropy= result_mi;
    cuts = 5;
    sep_mi = sep_misv(a_data,sorted_Spikes,Infos,seq ,cuts);
    analyse.sep_mi =  sep_mi;
    %add in spike trigger average options so it will be all in one file.
    
    
    %save analyse file
    cd(result_path)
    save( [analyse.filename(4:end) '_4'] ,'analyse')
    
    
end