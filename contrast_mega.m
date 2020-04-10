close all
clear all
clc

%this verison of spike are all zeroed
cd('D:\sfn2019\contrast\analyse\lump')
contrast = 1:5;
noise_level = 0;
rpt = 4;
mega = cell(1,length(contrast));
for con = 1:length(contrast)%contrast level
    
    sub=[];
    for rps = 1:rpt%repeats
        tcon = contrast(con);
        tnoi = 0;
        trps = rps;
        load(['c0' num2str(tcon) 'n' num2str(tnoi) '_' num2str(trps)])
        z = analyse;
        for i = 1:z.num_unit
            if z.staC(i)>300%this is the selection certia
                yo =length(sub)+1;
                sub(yo).sta_spikes = z.sta_spikes{1,i};
                sub(yo).sta_seq = z.sta_seq
                sub(yo).contrast = z.contrast;
                
                %sub(yo).type = typing(z.sta(:,i));
                sub(yo).sta = z.sta(:,i);
                sub(yo).spike_time = z.spike_time{1,i};
                sub(yo).spike_count = z.spike_count(i);
                sub(yo).id = z.id(i)*10+rps;
                sub(yo).ISI = z.ISI{1,i};
                sub(yo).TSMI = z.TSMI{1,i};
                sub(yo).t = z.mi_t;
                sub(yo).miph = z.miph(i);
                sub(yo).mipl = z.mipl(i);
                sub(yo).misv = z.misv(i);
                sub(yo).spkh = z.spike_entropy(i);
                sub(yo).sep_mi = z.sep_mi(:,i);
                sub(yo).seq= z.seq;
                sub(yo).condition = z.filename;
                analog_data  =  z.analog;
                SamplingRate=z.sampling_rate;
                norm_analog =analog_data-mean(analog_data);
                plot(norm_analog)
                pre_stamp = find(norm_analog>150,1);
                TimeStamp = (pre_stamp/SamplingRate) + 11;
                sub(yo).analog_stp= TimeStamp;
            end
        end
        
        
    end
    mega{1,con} = sub;
    
end