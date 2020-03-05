
%%% read file on a folder automaticlly 
name=dir('D:\Tina\tif*.grd')  %read all tif files under the folder
for ii=1:length(name);
    filename = char(name(ii,1).name);  %read the file name one by one
    c_filename = char(name(ii,1).name(1:8));  %read 1 to 8 char of the filename
    xyz = load(char(filename)); 
end

%%%read file with continus change file name
name = 'OG488_0309To0402_NB_' ;
for i = 1:1000
    filename = [name,num2str(i,'%02d'),'tif']
    %%% OG488_0309To0402_NB_0001
end