%%%% burst detect by Tina 2012, May
%%%% detect burst by spike interval < set_SpikeInterval and spike number > set_SpikeInBurst
%set_SpikeInterval = 0.25; % largest spike interval in a burst (s)
%set_SpikeInBurst = 5;  % smallist number of spikes in a burst

function [burst , BurstTime] = BurstDetect_Tina01(SpikeTime,set_SpikeInterval,set_SpikeInBurst)

%set_SpikeInterval = 0.25; % largest spike interval in a burst (s)
%set_SpikeInBurst = 5;  % smallist number of spikes in a burst


    aa = SpikeTime;
    b = length(aa) ;
    peak_int = aa(2:b)-aa(1:b-1) ; % interval of spikes
    a = peak_int < set_SpikeInterval ; % definition of burst : spike interval within 1 sec
    % note : 1111111101111111001111111011111  continous 1 means a burst, there are 4 bursts in this samples, at least 2 spikes can be a burst
    
    label_burst = bwconncomp(a) ; % use 'bwconncomp' to fine the number of burst
    label_list = label_burst.PixelIdxList ; % spikes index of every burst , type is 'cell'
    n_label =size(label_list) ; 
    n_label = n_label (2) ; % number of burst
    label_listOrg2 = label_list ;
    %%% remove the burst < n spikes
    for i = 1:n_label
        if length(label_list{i}) <set_SpikeInBurst 
        label_list{i} = [] ;
        end
    end
    c = cellfun('isempty',label_list);
    label_listOrg = label_list ;
    label_list(c>0.5) = [] ;

    n_label =size(label_list) ; 
    n_label = n_label (2);  % number of burst
    burst = {};
    BurstTime = [] ;
    
    for i = 1 : n_label
        bburst = aa(label_list{i}) ;
        burst(i) = {bburst} ; % spike time of every burst, type is 'cell'
%         BurstTime = [BurstTime aa(label_list{i})] ;
        BurstTime = horzcat(BurstTime ,bburst) ;
    end



