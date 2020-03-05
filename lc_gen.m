for i = 1:163
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

close all
lc = zeros(1,163);
comp2 = cell(2,6);
con = 1;
for i = 1:163
            figure;hold on
            for noi = 1:length(noise_level)%noise level
                z = mega{con,noi};
                plot(z(i).TSMI)             
            end
            hold off
      lc(i)= input('accept');
    close all
end
