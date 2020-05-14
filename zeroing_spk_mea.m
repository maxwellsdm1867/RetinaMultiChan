function opt =  zeroing_spk_mea(start_time,spikes)
for i = 1:size(spikes,2)
    opt{1,i}= zeroing_spk(start_time,spikes{1,i});
end
end
