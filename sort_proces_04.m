close all
clear all
clc

%load(uigetfile);% load file

%input condition and other parameter to
condition_name = ['542103524103'];%this is the noise level in order
cst_lev = ['222222333333'];%contrast level in order
num_trial = ['222222222222'];%trial number
name_seq = cell(1,length(condition_name));
for i = 1:length(condition_name)
    name_seq{i} = ['tr' num_trial(i) 'c0' cst_lev(i) 'n' condition_name(i)];
end

%loading sorted spikes
Spikes = cell(1,60);
channel = [12,13,14,15,16,17,21,22,23,24,25,26,27,28,31,32,33,34,35,36,37,38,41,42,43,44,45,46,47,48,51,52,53,54,55,56,57,58,61,62,63,64,65,66,67,68,71,72,73,74,75,76,77,78,82,83,84,85,86,87];
for h=1:60
    try
        Spikes{h} = eval(['adch_',int2str(channel(h))]);
    catch
        Spikes{h}=[];
    end
end

%count total units
total_units = 0;
for ck = 1:60
    if ~isempty(Spikes{1,ck})
        units=(Spikes{1,ck}(:,2))';
        total_units= total_units + length(unique(units));
    end
end

all_unit= {3,total_units};
s = 1;
for i = 1:60;
     if ~isempty(Spikes{1,i})
    time_stamp = (Spikes{1,i}(:,3))';
    unit_idx = (Spikes{1,i}(:,2))';
    num_subunit = length(unique(unit_idx));
    for j = 1:num_subunit
        temp_stamp  = time_stamp( unit_idx == j);
        if ~isempty(temp_stamp);
            all_unit{1,s}=temp_stamp ;
            all_unit{2,s} = j;
            all_unit{3,s} = channel(i);
            s= s+1;
        end
    end
     end
end


condition_time = 215.4;%seconds
trial_time = condition_time*ones(1,length(name_seq));
temp_unit = all_unit;
IDs = all_unit(2:3,:);
split_spike = cell(1,length(name_seq));
spacers = zeros(1,length(name_seq)+1);
for j = 1:length(name_seq)+1
    cnt = 0;
    for k = 1:j-1
        cnt = cnt+trial_time(k);
    end
    spacers(j) = cnt;
end

%ss = ['onoff' name_seq];
for i = 1:size(name_seq,2)
    vaild_range = [spacers(i) spacers(i+1)];
    tempC = cell(1,length(all_unit));
    for j = 1:length(all_unit)
        tk = temp_unit{1,j};
        tz =tk(tk>vaild_range(1));
        tempC{1,j}= tz(tz<=vaild_range(2))-vaild_range(1);
    end
    Spikes = [tempC; IDs];
    save(name_seq{i},'Spikes')
    
end

