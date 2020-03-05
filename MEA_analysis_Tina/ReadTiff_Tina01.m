AllImg = zeros(80,80,1000) ;
dir = 'D:\Tina\tif\'
name = 'OG488_0309To0402_NB_04' ;
for i = 1:1000
    filename = [dir,name,num2str(i,'%04d'),'.tif'] ;
    %%% OG488_0309To0402_NB_0001
    AllImg(:,:,i) = imread(filename) ;
end
SumIng = sum(AllImg,3) ; 
