% read mcd file 

%%%  check the streams type is ''Electrode Raw Data'



file = 'E:\Kevin\Retina analysis 20140417\20140411(MEA(3D5-21))frogretina1-6.mcd';

AllDataInfo =datastrm(file); % call for data info
t1=getfield(AllDataInfo,'sweepStartTime'); % call the start time
t2=getfield(AllDataInfo,'sweepStopTime') ; % call the stop time
StartStopVector = [t1 t2]; % time window
SamplingRate = getfield(AllDataInfo,'MillisamplesPerSecond2')/1000; % rescale sample rate to Hz
AllData = nextdata(AllDataInfo,'startend',StartStopVector,'originorder','on' ,'streamname','Electrode Raw Data'); % call whole data
rawdata = AllData.data(53,:); % call single channel
% timeline = [t1: 1/SamplingRate*1000:t2]/1000;
[vdata, tvals] = ad2muvolt(AllDataInfo, rawdata,  getfield(AllDataInfo,'StreamNames')); % rescale vdata(uV) , tvals(us)
tvals = tvals/1E6; % rescale unit to s
plot(tvals,vdata);

data = vdata ;
[b,a] = butter(2,200/10000,'high'); % set butter filter
FilterData = filter(b,a,vdata);
figure(2),subplot(211),plot(tvals,vdata),subplot(212),plot(tvals,FilterData)
std_data = std(FilterData);
[pks,locs] = findpeaks(FilterData,'MINPEAKHEIGHT', std_data*5);
plot(-FilterData), hold on, plot(locs,std_data*5,'o'),hold off 