classdef SR_analysis
    %This is the class for SR lump data analysis
    %   After spike sorting and basic analysis, this class is a way of
    %   collecting the functions and the datas together that is related to
    %   the statical analyis
    
    properties
        mega
        lc
        
    end
    properties  (Access = public)
        %Parameters that we feed into the function handle
        fr_BinningInterval = 0.040;%firing rate intercval
        mi_BinningInterval = 0.040;%time shift MI window
        num_cuts = 10;%number of cuts in the cut_handel(obj)
        i_cut = 1;%ith cut in the cut_handel(obj)
        i_side = 1;%ith slide in sliding_window(obj)
        slide_width = 100;%width of one slide in sliding_window(obj)
        win_width = 3000;%width of a window in sliding_window(obj)
        max_range = [140 160];
        DataTime = 200;%data time which is fixed
        num_unit = 271;
        cell_idx = [];
    end
    properties  (Access = private)
        %Parameters that come with the experiment
        stim_rate = 5;%update rate of the stimulus
        led_rate = 25;%update rate of the led
        num_contrast = 2;%total contrast level
        num_noise = 6;%total number of the noise level
    end
    methods
        function obj = SR_analysis(mega,lc)
            % Construct an instance of this class
            %   if the number of the input is 2 then use this constructor
            if nargin == 2
                obj.mega = mega;
                obj.lc = lc;
            end
        end
        
        
        %% extractors, to extract basic information form mega file according to
        function rslt= sr_extracter_sv(obj,varable_nme,cell_tpe)
            %This is the function to extract single value variable cell by cell
            %  the result is an 1 by 2 cell and each cell contains the
            cdt_c = size(obj.mega,1);
            cdt_n = size(obj.mega,2);
            rslt = cell(cdt_c,1) ;
            for i = 1:obj.num_unit
                for con = 1:cdt_c%contrast level
                    temp = [];
                    for noi = 1:cdt_n%noise level
                        z = obj.mega{con,noi};
                        if obj.lc(i) == 1&& ismember(z(i).type,cell_tpe)
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
            for i = 1:obj.num_unit
                for con = 1:cdt_c%contrast level
                    temp = {};
                    for noi = 1:cdt_n%noise level
                        z = obj.mega{con,noi};
                        
                        if obj.lc(i) == 1&& ismember(z(i).type,cell_tpe)%ismember(obj.cell_idx(i),cell_tpe)
                            
                            ks1 =  eval(['z(' num2str(i) ').' varable_nme]);
                            temp{size(temp,1)+1,1}= ks1;
                        end
                    end
                    %keyboard
                    rslt{con,1} = [rslt{con,1} temp];
                end
            end
        end
        %% post clustering extraction
        function rslt= sr_extracter_sv_p(obj,varable_nme,cell_tpe)
            %This is the function to extract single value variable cell by cell
            %  the result is an 1 by 2 cell and each cell contains the
            cdt_c = size(obj.mega,1);
            cdt_n = size(obj.mega,2);
            rslt = cell(cdt_c,1) ;
            for i = 1:obj.num_unit
                for con = 1:cdt_c%contrast level
                    temp = [];
                    for noi = 1:cdt_n%noise level
                        z = obj.mega{con,noi};
                        if obj.lc(i) == 1&&  ismember( obj.cell_idx(i),cell_tpe)
                            ks1 =  eval(['z(' num2str(i) ').' varable_nme]);
                            temp = [temp; ks1];
                        end
                    end
                    rslt{con,1} = [rslt{con,1} temp];
                end
            end
            
        end
        
        function rslt= sr_extracter_mv_p(obj, varable_nme,cell_tpe)
            %Extract the non-single value data form mega
            %   input is the object it self and what variable you want to
            %   extract, and what specific cell type you want extract
            
            cdt_c = size(obj.mega,1);
            cdt_n = size(obj.mega,2);
            rslt = cell(cdt_c,1) ;
            for i = 1:obj.num_unit
                for con = 1:cdt_c%contrast level
                    temp = {};
                    for noi = 1:cdt_n%noise level
                        z = obj.mega{con,noi};
                        
                        if obj.lc(i) == 1&& ismember(obj.cell_idx(i),cell_tpe)
                            
                            ks1 =  eval(['z(' num2str(i) ').' varable_nme]);
                            temp{size(temp,1)+1,1}= ks1;
                        end
                    end
                    %keyboard
                    rslt{con,1} = [rslt{con,1} temp];
                end
            end
        end
        %% function mapper  of different type
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
        %%
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
            fr = @(x) BinSpk1(obj.fr_BinningInterval,x,obj.DataTime);
            function [BinningSpike] = BinSpk1(BinningInterval,spike_time,DataTime)
                % transfer spike time into firing rate
                %   input binning interval in seconds, spike in seconds
                BinningTime = [0 :  BinningInterval : DataTime];
                %[n,xout] = hist(spike_time,BinningTime);  %putting spikes in the right timings according to an assigned bin interval
                [N,edges] = histcounts(spike_time,BinningTime);%updated verision of binning spikes
                BinningSpike = N;
            end
            
        end
        
        function cut_mi = cut_mi_handel(obj)
            cut_mi = @(x,y) cut_mis(x,y,obj.num_cuts);
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
        
        
        function cut = cut_handel(obj)
            %cut the array in this class;
            cut = @(p) spiltter(obj.num_cuts,p,obj.i_cut);
            function rslt = spiltter(num_cuts,n,i_cuts)
                ed = (length(n)/num_cuts)*i_cuts;
                st = ed-(length(n)/num_cuts)+1;
                rslt = n(st:ed);
            end
            
        end
        
        function res = sliding_window_handel(obj)
            %create a sliding window
            res = @(p)  slider(p,obj.slide_width,obj.win_width, obj.i_side);
            function rslt = slider(x,slide_width,win_width,i_win)
                rslt =x(((i_win-1)*slide_width+1):(((i_win-1)*slide_width)+win_width-1));
            end
        end
        
        function mi = tsmi_handle(obj)
            mi = @(x,y) find_MI(x,y,(1/obj.mi_BinningInterval));
            function rslt = find_MI(sys_opt,seq,BinningSamplingRate)
                %[MI, MI_shuffled,t] = tsmi_clean1(sys_opt,seq);
                [MI,t]= only_timeshift(seq,sys_opt,BinningSamplingRate);
                rslt = MI;
            end
        end
        
        function rslt = max_handle(obj,x)
            rslt =@(x) max(x(min(obj. max_range):max(obj. max_range)));
        end
        
        function rslt = mean_handle(obj,x)
            rslt = @(x) mean(x);
        end
        
        function rslt = err_handle(obj,x)
            rslt =  @(x) std(x);%/length(x);
        end
        
        function vis_all_mv(obj,dat)
            %multi varible vaiulization
            num_cell = size(dat{1,1},2);
            cdt = 6;
            for i = 1:num_cell
                figure
                subplot(2,1,1)
                hold on
                for j = 1:cdt
                    plot(dat{1,1}{j,i})
                end
                hold off
                subplot(2,1,2)
                hold on
                for k = 1:cdt
                    plot(dat{2,1}{k,i})
                end
                hold off
            end
        end
        
        function rc = lc_gen(obj,dat)
            %for generate the lc, after the mega is generated
            num_cell = size(dat{1,1},2);
            rc = zeros(1,num_cell);
            cdt = 6;
            for i = 1:num_cell
                figure
                subplot(2,1,1)
                hold on
                for j = 1:cdt
                    plot(dat{1,1}{j,i})
                end
                hold off
                subplot(2,1,2)
                hold on
                for k = 1:cdt
                    plot(dat{2,1}{k,i})
                end
                hold off
                rc(i)= input('accept');
                close all
            end
            
        end
        function rslt = collapse(obj,dat)
            %collapse the data of either silding window or cut to a single
            %cell, to use other builting in this class
            num_rpt = length(dat);%number of window or cuts
            num_cell = size(dat{1,1}{1,1},2);%number of cell in the specified type
            c = obj.num_contrast;%number of noise
            n = obj.num_noise;
            temp = cell(2,1);
            
            for i = 1:num_cell%cell number
                
                for j = 1:c%contrast levels
                    for k = 1:n%noise levels
                        holder = zeros(1,10);
                        for m = 1: num_rpt%number repeats
                            holder(1,m) =  dat{1,m}{j,1}{k,i};
                        end
                        temp{j,1}{k,i} = holder;
                    end
                end
                
            end
            rslt = temp;
        end
        
        %% code for linear nonliner extraction from gaussain stimulus
        function sta = sta_handle(obj)
            sta= @(x,y) find_sta(x,y);
            function rslt = find_sta(stim,fr)
                %[MI, MI_shuffled,t] = tsmi_clean1(sys_opt,seq);
                [stas,~] = obj.sta_check(stim,fr);
                rslt = stas;
            end
        end
        function sta =gen_sg_handle(obj)
            sta= @(x,y) find_sta(x,y);
            function rslt = find_sta(stim,fr)
                %[MI, MI_shuffled,t] = tsmi_clean1(sys_opt,seq);
                [~,gen_signal,~,~] = obj.sta_check(stim,fr);
                rslt = gen_signal;
            end
        end
        
        function sta =nlx_handle(obj)
            sta= @(x,y) find_sta(x,y);
            function rslt = find_sta(stim,fr)
                %[MI, MI_shuffled,t] = tsmi_clean1(sys_opt,seq);
                [~,~,nlx,~] = obj.sta_check(stim,fr);
                rslt = nlx;
            end
        end
        
        function sta =nly_handle(obj)
            sta= @(x,y) find_sta(x,y);
            function rslt = find_sta(stim,fr)
                %[MI, MI_shuffled,t] = tsmi_clean1(sys_opt,seq);
                [~,~,~,nly] = obj.sta_check(stim,fr);
                rslt = nly;
            end
        end
        
        function [sta,gen_signal,fx,fy] = sta_check(~,sta_seq,sta_fr)
            sps = sta_fr;
            nsp = sum(sps);
            seq1 = sta_seq(1:7500);
            Stim= seq1-mean(seq1);
            Stim = Stim';
            sps = sps';
            ntfilt = 13;  % length of linear the filter
            paddedStim = [zeros(ntfilt-1,1); Stim]; % pad early bins of stimulus with zero
            Xdsgn = hankel(paddedStim(1:end-ntfilt+1), Stim(end-ntfilt+1:end));%hankel matrix
            sta = (Xdsgn'*sps)/nsp;
            %gen_signal =(Xdsgn)*sta; %the output of sta(Generator output)
            rawfilteroutput0 = conv(Stim,flip(sta)) ;
            gen_signal  = rawfilteroutput0(1:7500);
            %scatter(gen_signal,sps) make a scatter handle later
            nfbins = 11;
            [cts,binedges,binID] = histcounts(gen_signal,nfbins);
            fx = binedges(1:end-1)+diff(binedges(1:2))/2; % use bin centers for x positions
            fy = zeros(nfbins,1); % y values for nonlinearity
            for jj = 1:nfbins
                fy(jj) = mean(sps(binID==jj));
            end
            
        end
        
        %% code for cell type clustering
        function sta =res_handle(obj)%biphasicness of the filter
            sta= @(x,y) find_celltype(x,y);
            function rslt = find_celltype(sta,nly)
                %[MI, MI_shuffled,t] = tsmi_clean1(sys_opt,seq);
                [res,u_idx,r]= cell_type_pram(obj,sta,nly);
                rslt = res;
            end
        end
        function sta =u_idx_handle(obj)% u-shanpe of the nonlinearity
            sta= @(x,y) find_celltype(x,y);
            function rslt = find_celltype(sta,nly)
                %[MI, MI_shuffled,t] = tsmi_clean1(sys_opt,seq);
                [res,u_idx,r]= cell_type_pram(obj,sta,nly);
                rslt = u_idx;
            end
        end
        function sta =r_handle(obj)%ploarity
            sta= @(x,y) find_celltype(x,y);
            function rslt = find_celltype(sta,nly)
                %[MI, MI_shuffled,t] = tsmi_clean1(sys_opt,seq);
                [res,u_idx,r]= cell_type_pram(obj,sta,nly);
                rslt = r;
            end
        end
        function [res,u_idx,r]= cell_type_pram(obj,sta,nly)
            fact = 100;
            up_sta= interp(sta,fact);
            f_up= up_sta(1:1200);
            u_idx = stm_ushape(nly);
            f_up= up_sta(1:1200);
            temp = f_up(800:end);
            m = mean(f_up(1:600));
            mid = temp-m;
            po = abs(sum(mid(mid>0)));
            ne = abs(sum(mid(mid<0)));
            tot_area = sum(abs(mid));
            d= po-ne;
            res = d./tot_area;
            zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
            pr= zci(f_up);
            z = [];
           
            %need to hadle the points without zeros crossing
            if sum(abs(f_up))==abs(sum(f_up))
                r = sum(f_up);
            else%one with zero crossing
                for i= 1:length(pr)
                    temp = pr(i);
                    z = [z sum(f_up(temp:end))];
                end
                [m,mi] = max(abs(z));
                r = z(mi);
            end
            %             ver 2
            %             fact = 100;
            %             up_sta= interp(sta,fact);
            %             f_up= up_sta(1:1200);
            %             u_idx = stm_ushape(nly);
            %             f_up= up_sta(1:1200);
            %             temp = f_up(800:end);
            %             m = mean(f_up(1:600));
            %             mid = temp-m;
            %             po = abs(sum(mid(mid>0)));
            %             ne = abs(sum(mid(mid<0)));
            %             tot_area = sum(abs(mid));
            %             d= po-ne;
            %             res = d./tot_area;
            %
            %             r = sum(f_up(1000:end));
        end
        %% others
        function [ps,hs] = sr_condtion(obj,dat,m)
            %m is the mean matrix
            %dat is the colapse matrix
            %this function only condiers weather noise can increse the tsmi
            num_cell = size(dat{1,1},2);%number of cell in the specified type
            ps = zeros(2,num_cell);
            hs = zeros(2,num_cell);
            for i = 1:obj.num_contrast
                for j = 1:num_cell
                    
                    temp_a = m{i,1}(:,j);%all of them for index purpose
                    [v,po] = max( temp_a);
                    po
                    %max value is not one and it is significent after the
                    %sign rank test.
                    if po~=1
                        zn = dat{i,1}{1,j};%zero noise conditions
                        comp =  dat{i,1}{po,j};
                        [p,h] = signrank(zn,comp);
                        ps(i,j) = p;
                        hs(i,j) = h;
                    end
                end
            end
        end
        
        function [rslt]= row_extracter(obj,data,con,noi)
            %con is the constrast level and noi is the noise level
            temp =data{con,1};
            rslt = temp(noi,:);
        end
        
    end%methods ends here
    
    
    
    
    
    
    
    
    
end