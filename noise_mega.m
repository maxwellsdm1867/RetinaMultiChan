close all
clear all
clc
cd('D:\sf\noise_lump')
contrast = 2:3;
noise_level = 0:5;
rpt = 4;
mega = cell(length(contrast),length(noise_level));
for con = 1:length(contrast)%contrast level
    for noi = 1:length(noise_level)%noise level
        sub=[];
        for rps = 1:rpt%repeats
            tcon = contrast(con);
            tnoi = noise_level(noi);
            trps = rps;
            load(['c0' num2str(tcon) 'n' num2str(tnoi) '_' num2str(trps)])
            z = analyse;
            for i = 1:z.num_unit
                if z.staC(i)>300%this is the selection certia
                    yo =length(sub)+1;
                    sub(yo).type = typing(z.sta(:,i));
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
                    sub(yo).condition = z.filename(4:end);
                     sub(yo).Stim = z.Stim;
                    sub(yo).staSpikes = z.staSpikes{1,i};
                    sub(yo).staSeq = z.staSeq;
                    sub(yo).nlx = z.nlx{1,i};
                    sub(yo).nly = z.nly{1,i};
                    sub(yo).nlfit = z.nlfit{1,i};
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
        mega{con,noi} = sub;
    end
end