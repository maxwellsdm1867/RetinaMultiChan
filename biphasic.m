%test the bipashic index, first we upasample the sta for 7th to 13th
pre_up = pre_clust;
%upsample by a factor of 200
fact = 100;
post_up = zeros(size(pre_up,1)*fact,size(pre_up,2));
for i = 1:size(pre_up,2)
    post_up(:,i) = interp(pre_up(:,i),fact);
end
f_up= post_up(1:1200,:);
figure; hold on
for i=1:size(pre_clust,2)
    
    plot(f_up(:,i))
    %plot(pre_clust(:,i))
    %title(num2str(test(i)))
    
end
hold off

temp = f_up(800:end,:);
m = repmat(mean(f_up(1:600,:)),size(temp,1),1);
mid= temp-m;
po = abs(sum(mid.*(mid>0)));
ne = abs(sum(mid.*(mid<0)));
tot_area = sum(abs(mid));
d= po-ne;
res = d./tot_area;
figure; hold on
for i=1:size(pre_clust,2)
    if res(i)<0.8 & res(i)>-0.8 & r(i)==1
        plot(f_up(:,i))
        ylim([0.027 0.038])
    end
end
hold off

%onoff or offon use the derivative
com = mid(100:end,:)<0;
k = diff(com);
r = zeros(1,size(k,2));
for i = 1:size(k,2)
    temo = k(:,i);
    aa =find(temo~=0);
    if~isempty(aa)
        r(i) = temo(aa(1));
    end
end
