for i = 1:92
    for con = 1:length(contrast)%contrast level
        for noi = 1:length(noise_level)%noise level
            z = mega{con,noi};
            if z(i).type == tar 
                tempm(con,noi) = mean(z(i).sep_mi);
                temps(con,noi) = std(z(i).sep_mi)/sqrt(5);
                tar_idx(i)=1;
                comp1{4,i}=z(i).id;
            else
                tar_idx(i) = 0;
                comp1{3,i}=0;
            end
        end
    end
end

%%
close all 
num_unit = 163;
noise_level = 0:5;
lc = zeros(1,num_unit);
comp2 = cell(2,6);
con =2;
for i = 1:num_unit
            figure;hold on
            for noi = 1:length(noise_level)%noise level
                z = mega{con,noi};
                plot(z(i).t,z(i).TSMI)
                xlim([-2000 2000])
            end
            hold off
      lc(i)= input('accept');
    close all
end
