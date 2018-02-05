% Test the saturation detection
close all
clear 
clc

demdir=['./'];
%fprintf('working: %s\n',demdir);
%f = dir([demdir,'/*.tif']);

f{1}='ex1.tif';
res='8'; % Note the data type must be character, with value less than 16. '0.517'

display(['Number of files: ',num2str(length(f)),''])
dataratio=[];
for i=1:length(f)  
    close all
 
    matchFile = [demdir,'/',f{i}];
    fprintf('processing %d of %d: %s \n',i,length(f),matchFile)
    
    infile= matchFile;
    
%   data=readGeotiff(infile);
   data.z=imread(infile);
    
    [ms,ratio] = SatuFilter(data.z,res,1); % takes 3 minutes for 8m images
    % Saturation ratio defined as the total number of saturated points
    % divided by the total number of data points
    if ratio > 1
       display(['This image is saturated, with saturation ratio as ',num2str(ratio),'%.'])
    else
       display(['This image is not saturated.'])
    end
   
    dataratio=[dataratio;{num2str(i),num2str(ratio),infile}];
    otype=['./'];
    
    h=double(data.z);%h=I.z;h=m.z;
    minh=min(min(h(h~=0)));maxh=max(max(h(h~=0)));
    meanh=mean(mean(h(h~=0)));stdh=std(h(h~=0));
    figure; imagesc(data.z);caxis([minh maxh*1.3]);colorbar;colormap jet; %chunli dai
    title(['Input image']);
    caxis([max(meanh-stdh*0.5, minh) min(meanh+1*stdh,maxh)]);colorbar;%
    OutOrthoName= strrep(f{i},'.tif','_browse.tif');
    OutOrthoName=[otype,OutOrthoName];
    print('-dtiff','-r200',OutOrthoName)
    
    OutSatuMaskName= strrep(f{i},'.tif','_satumask.tif');
    OutSatuMaskName=[otype,OutSatuMaskName];   
    figure; imagesc(double(ms));colormap gray;colorbar;caxis([0 1]) %chunli dai
    title(['Saturation mask']);
    print('-dtiff','-r500',OutSatuMaskName)

end
% dataratio
[nrows,ncols] = size(dataratio);
fid = fopen('dataratio.txt', 'a');
for row = 1:nrows
    fprintf(fid,'%s %s %s\n',dataratio{row,:});
end
fclose(fid);

