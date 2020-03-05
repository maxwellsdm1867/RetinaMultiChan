classdef noise_analysis
    %This is the class for SR lump data analysis
    %   After spike sorting and basic analysis, this class is a way of
    %   collecting the functions and the datas together that is related to
    %   the statical analyis
    
    properties
        mega
        lc
        
    end
    properties  (Access = private)
        BinningInterval = 0.040;
        DataTime = 200;
        
    end
    methods
        function obj = noise_analysis(mega,lc)
            % Construct an instance of this class
            %   if the number of the input is 2 then use this constructor
            if nargin == 2
                obj.mega = mega;
                obj.lc = lc;
            end
        end
        
        function rslt= sr_extracter_sv(obj,varable_nme,cell_tpe)
            %This is the function to extract single value variable cell by cell
            %  the result is an 1 by 2 cell and each cell contains the
            cdt_c = size(obj.mega,1);
            cdt_n = size(obj.mega,2);
            rslt = cell(cdt_c,1) ;
            for i = 1:163
                for con = 1:cdt_c%contrast level
                    temp = [];
                    for noi = 1:cdt_n%noise level
                        z = obj.mega{con,noi};
                        if obj.lc(i) == 1&& z(i).type == cell_tpe
                            ks1 =  eval(['z(' num2str(i) ').' varable_nme]);
                            temp = [temp; ks1];
                        end
                    end
                    rslt{con,1} = [rslt{con,1} temp];
                end
            end
            
        end
        
        function rslt= sr_extracter_mv(obj, varable_nme,cell_tpe)
            %Extract the non-single value data form mega
            %   input is the object it self and what variable you want to
            %   extract, and what specific cell type you want extract
            
            cdt_c = size(obj.mega,1);
            cdt_n = size(obj.mega,2);
            rslt = cell(cdt_c,1) ;
            for i = 1:163
                for con = 1:cdt_c%contrast level
                    temp = {};
                    for noi = 1:cdt_n%noise level
                        z = obj.mega{con,noi};
                        if obj.lc(i) == 1&& z(i).type == cell_tpe
                            ks1 =  eval(['z(' num2str(i) ').' varable_nme]);
                            temp{size(temp,1)+1,1}= ks1;
                        end
                    end
                    %keyboard
                    rslt{con,1} = [rslt{con,1} temp];
                end
            end
        end
        
        
        
        % these are function appliers that is either singular or binary
        %singular means it takes one input and map to one output using the given figure handle
        function rslt = singular_mv2sv(obj,input_cel,function_handle)
            %this function will map the function that takes single value of the data
            %and apply it to every posotion of the array, those function all have
            %customized functions that only take single values
            
            
            rslt = cell(size(input_cel,1),size(input_cel,2));
            for s = 1:size(input_cel,1)
                temp = input_cel{s,1};
                temp_rslt = zeros(size(temp,1),size(temp,2));
                for i = 1:size(temp,1)
                    for j = 1:size(temp,2)
                        temp_rslt(i,j)= function_handle(temp{i,j});
                    end
                end
                rslt{s,1}= temp_rslt;
            end
            
            
        end
        function rslt = singular_mv2mv(obj,input_cel,function_handle)
            %this function will map the function that takes single value of the data
            %and apply it to every posotion of the array, those function all have
            %customized functions that only take single values
            
            
            rslt = cell(size(input_cel,1),size(input_cel,2));
            for s = 1:size(input_cel,1)
                temp = input_cel{s,1};
                temp_rslt = cell(size(temp,1),size(temp,2));
                for i = 1:size(temp,1)
                    for j = 1:size(temp,2)
                        temp_rslt{i,j}= function_handle(temp{i,j});
                    end
                end
                rslt{s,1}= temp_rslt;
            end
            
            
        end
        function rslt = binary_mv2sv(obj,input_cel1,input_cel2,function_handle)
            %this function will map the given non-single value function to
            
            
            rslt = cell(size(input_cel1,1),size(input_cel1,2));
            for s = 1:size(input_cel1,1)
                temp1 = input_cel1{s,1};
                temp2 = input_cel2{s,1};
                temp_rslt = zeros(size(temp1,1),size(temp1,2));
                for i = 1:size(temp1,1)
                    for j = 1:size(temp1,2)
                        temp_rslt(i,j)= function_handle(temp1{i,j},temp2{i,j});
                    end
                end
                rslt{s,1}= temp_rslt;
            end
            
            
        end
        function rslt = binary_mv2mv(obj,input_cel1,input_cel2,function_handle)
            
            rslt = cell(size(input_cel1,1),size(input_cel1,2));
            for s = 1:size(input_cel1,1)
                temp1 = input_cel1{s,1};
                temp2 = input_cel2{s,1};
                temp_rslt = cell(size(temp1,1),size(temp1,2));
                for i = 1:size(temp1,1)
                    for j = 1:size(temp1,2)
                        temp_rslt{i,j}= function_handle(temp1{i,j},temp2{i,j});
                    end
                end
                rslt{s,1}= temp_rslt;
            end
            
            
        end
        function rslt = binary_mix2sv(obj,input_sv,input_mv,function_handle)
            
            rslt = cell(size(input_sv,1),size(input_sv,2));
            for s = 1:size(input_sv,1)
                temp_sv = input_sv{s,1};
                temp_mv = input_mv{s,1};
                temp_rslt = zeros(size(temp1,1),size(temp1,2));
                for i = 1:size(temp1,1)
                    for j = 1:size(temp1,2)
                        temp_rslt(i,j)= function_handle(temp_sv(i,j),temp_mv{i,j});
                    end
                end
                rslt{s,1}= temp_rslt;
            end
            
            
        end
        function rslt = binary_mix2mv(obj,input_sv,input_mv,function_handle)
            
            rslt = cell(size(input_sv,1),size(input_sv,2));
            for s = 1:size(input_sv,1)
                temp_sv = input_sv{s,1};
                temp_mv = input_mv{s,1};
                temp_rslt = cell(size(temp_sv,1),size(temp_sv,2));
                for i = 1:size(temp_sv,1)
                    for j = 1:size(temp_sv,2)
                        temp_rslt{i,j}= function_handle(temp_sv(i,j),temp_mv{i,j});
                    end
                end
                rslt{s,1}= temp_rslt;
            end
            
            
        end
        
        %these function bloew will be functions that create the modified function that can
        %be work with sr_applier_mv
        % for a mixed binary, single value is the first input and
        % non-single valued is the latter, these are all handel generators,
        % which includes their own functions
        function adj_spk = adj_spk_hadel(obj)
            %binary_mv2mv
            adj_spk = @(sv,mv) zeroing_spk(sv,mv);
            function opt =  zeroing_spk(start_time,spikes)
                TimeStamp = start_time;
                spike_time = spikes;
                for n = 1:length(spike_time)
                    if (spike_time(n)- TimeStamp) < 0
                        spike_time(n) = -1;
                    else
                        spike_time(n) =  spike_time(n)-TimeStamp;
                    end
                end
                opt = spike_time(spike_time>=0);
            end
        end
        function  fr = fr_handel(obj)
            %singular_mv2mv
            %single value handel for the binning spike function
            fr = @(x) BinSpk1(obj.BinningInterval,x,obj.DataTime);
            function [BinningSpike] = BinSpk1(BinningInterval,spike_time,DataTime)
                % transfer spike time into firing rate
                %   input binning interval in seconds, spike in seconds
                BinningTime = [0 :  BinningInterval : DataTime-BinningInterval];
                [n,xout] = hist(spike_time,BinningTime);  %putting spikes in the right timings according to an assigned bin interval
                BinningSpike = n;
                BinningSpike(:,1) = 0;BinningSpike(:,end) = 0;
            end
            
        end
        function tsmi = tsmi_handel(obj)
            %mv2mv binary
            tsmi = @(resp,stim) tsmi_clean(resp,stim);
            function MI = tsmi_clean(sys_opt,seq)
                %For computating the TSMI curve between any 2 equal length sequence
                %   It requires two sequence to have the samle size(sampling rate)
                states = 12;
                shift = [6000 6000];
                BinningSamplingRate = 40;
                BinningInterval = 1/BinningSamplingRate;
                bin = BinningInterval*1000;
                X =seq;
                nX = sort(X);
                abin = length(nX)/states;
                intervals = [ nX(1:abin:end) inf];  %make the states equally distributed
                for k = 1:length(X)
                    [a,b] = find(X(k)<intervals,1);
                    isi2(k) = b-1;  %a new inter-pulse-interval denoted by the assigned states (ex:1-5) and with the same sampling rate as the BinningSpike
                end
                
                Neurons = sys_opt;%response
                %Neurons=isi2;
                backward=ceil(shift(1)/bin); forward=ceil(shift(2)/bin);
                dat=[];informationp=[];temp=backward+2;
                for i=1:backward+1 %past(t<0)
                    x = Neurons((i-1)+forward+1:length(Neurons)-backward+(i-1))';
                    y = isi2(forward+1:length(isi2)-backward)';
                    dat{i}=[x,y];
                    [N,C]=hist3(dat{i}); %20:dividing firing rate  6:# of stim
                    px=sum(N,1)/sum(sum(N)); % x:stim
                    py=sum(N,2)/sum(sum(N)); % y:word
                    pxy=N/sum(sum(N));
                    temp2=[];
                    for j=1:length(px)
                        for k=1:length(py)
                            temp2(k,j)=pxy(k,j)*log( pxy(k,j)/ (py(k)*px(j)) )/log(2);
                        end
                    end
                    temp=temp-1;
                    informationp(temp)=nansum(temp2(:));
                end
                
                dat=[];informationf=[];temp=0;sdat=[];
                
                for i=1:forward
                    x = Neurons(forward+1-i:length(Neurons)-backward-i)';
                    y = isi2(forward+1:length(isi2)-backward)';
                    dat{i}=[x,y];
                    
                    [N,C]=hist3(dat{i}); %20:dividing firing rate  6:# of stim
                    px=sum(N,1)/sum(sum(N)); % x:stim
                    py=sum(N,2)/sum(sum(N)); % y:word
                    pxy=N/sum(sum(N));
                    temp2=[];
                    for j=1:length(px)
                        for k=1:length(py)
                            temp2(k,j)=pxy(k,j)*log( pxy(k,j)/ (py(k)*px(j)) )/log(2);
                        end
                    end
                    temp=temp+1;
                    informationf(temp)=nansum(temp2(:));
                end
                information=[informationp informationf];
                t=[-backward*bin:bin:forward*bin];
                
                
                MI = information(151);
                
%                 % shuffle
%                 r=randperm(length(Neurons));
%                 for j=1:length(r)
%                     sNeurons(j)=Neurons(r(j));
%                 end
%                 Neurons=sNeurons;
%                 backward=ceil(shift(1)/bin); forward=ceil(shift(2)/bin);
%                 dat=[];informationp=[];temp=backward+2;
%                 for i=1:backward+1 %past(t<0)
%                     x = Neurons((i-1)+forward+1:length(Neurons)-backward+(i-1))';
%                     y = isi2(forward+1:length(isi2)-backward)';
%                     dat{i}=[x,y];
%                     [N,C]=hist3(dat{i}); %20:dividing firing rate  6:# of stim
%                     px=sum(N,1)/sum(sum(N)); % x:stim
%                     py=sum(N,2)/sum(sum(N)); % y:word
%                     pxy=N/sum(sum(N));
%                     temp2=[];
%                     for j=1:length(px)
%                         for k=1:length(py)
%                             temp2(k,j)=pxy(k,j)*log( pxy(k,j)/ (py(k)*px(j)) )/log(2)/(bin/1000);
%                         end
%                     end
%                     temp=temp-1;
%                     informationp(temp)=nansum(temp2(:));
%                 end
%                 
%                 dat=[];informationf=[];temp=0;sdat=[];
%                 
%                 for i=1:forward
%                     x = Neurons(forward+1-i:length(Neurons)-backward-i)';
%                     y = isi2(forward+1:length(isi2)-backward)';
%                     dat{i}=[x,y];
%                     
%                     [N,C]=hist3(dat{i}); %20:dividing firing rate  6:# of stim
%                     px=sum(N,1)/sum(sum(N)); % x:stim
%                     py=sum(N,2)/sum(sum(N)); % y:word
%                     pxy=N/sum(sum(N));
%                     temp2=[];
%                     for j=1:length(px)
%                         for k=1:length(py)
%                             temp2(k,j)=pxy(k,j)*log( pxy(k,j)/ (py(k)*px(j)) )/log(2)/(bin/1000);
%                         end
%                     end
%                     temp=temp+1;
%                     informationf(temp)=nansum(temp2(:));
%                 end
%                 information=[informationp informationf];
%                 MI_shuffled = information;%unshuffle minous shuffled
                
            end
            
            
        end
        function cut_mi = cut_mi_handel(obj)
            cut_mi = @(x,y) cut_mis(x,y,10);
            function rslt = cut_mis(seq,BinningSpike,cuts)
                
                isi = seq;
                rslt =[];
                states =16;
                %% Stimuli process
                for i = 1:cuts
                    Neurons = BinningSpike(1,1+((i-1)*length(BinningSpike)/cuts):i*length(BinningSpike)/cuts);
                    X =isi;%(1,1+((i-1)*length(BinningSpike)/cuts):i*length(BinningSpike)/cuts);
                    nX = sort(X);
                    abin = length(nX)/states;
                    intervals = [ nX(1:abin:end) inf];  %make the states equally distributed
                    for k = 1:length(X)
                        [a,b] = find(X(k)<intervals,1);
                        isi2(k) = b-1;  %a new inter-pulse-interval denoted by the assigned states (ex:1-5) and with the same sampling rate as the BinningSpike
                    end
                    x = Neurons;
                   % y = isi2;
                    y = isi2(1,1+((i-1)*length(BinningSpike)/cuts):i*length(BinningSpike)/cuts);%post spilcing
                    dat=[x;y];
                    [N,C]=hist3(dat'); %20:dividing firing rate  6:# of stim
                    px=sum(N,1)/sum(sum(N)); % x:stim
                    py=sum(N,2)/sum(sum(N)); % y:word
                    pxy=N/sum(sum(N));
                    lpx = -log2(px);
                    lpy = -log2(py);
                    lpxy = -log2(pxy);
                    temp = nansum(px.*lpx)+nansum(py.*lpy)-sum(nansum(pxy.*lpxy));
                    rslt = [rslt temp];
                end
            end
        end
        function cplot_1(dat,er)
            keyborad
            num_cell = size(dat{1,1},2);
            for i = 1:num_cell
                figure
                subplot(1,2,1)
                errorbar(dat{1,1}(:,i),er{1,1}(:,i))
                subplot(1,2,2)
                errorbar(dat{2,1}(:,i),er{2,1}(:,i))
            end
        end
        
        
        
        
    end
    
    
    
    
    
    
    
    
    
    
    
end