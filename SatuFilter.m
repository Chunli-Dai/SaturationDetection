function [M,ratio]= SatuFilter(varargin)
% [M,ratio]= SatuFilter(m0,res,opt) 
% SatuFilter returns a mask for the input image using wavelet filter.
%  The mask identifies the location of saturation by a wavelet transform,
%  used to detect periodic variations of brightness (i.e. striping) with
%  varying frequencies at different locations, in raw, orthorectified imagery.
%
%  Output variables:
%    M:    the 2D saturation mask, 1 is saturated, 0 is not saturated.
%    ratio: the saturation ratio, which is the percentage of points that 
% 	    are saturated over the total points.
%  Input variables:
%    m0:  the 2D matrix of raw image or orthoimage.
%    res: the resolution of input image, in character array, e.g. '8'
%    opt: a control parameter for options
%         opt=1 for raw images, and opt=0 for orthoimages
%         opt=1, means the direction of stripes is known as 90 degree.
%	      otherwise, the direction needs to be searched.
%
% Chunli Dai, dai.56@buckeyemail.osu.edu
% Version 1. April, 2017
%

%Initialize input variables
m0 = varargin{1};
resc= varargin{2};res=str2double(resc);
opt=varargin{3};

cksatu1=clock;

%Set control parameters
resr=8;%16;% Reduced resolution in meter;
% be careful, this resolution cannot be too rough. Suggested value 8, 16, 24, 32
% main parameters
Ts1=680./(resr); % % the mean of stripes wavelength
Ts1min=280./(resr);Ts1max=1600./(resr); %  min and maximum of the stripes wavelength
dsp=round(240./(resr)); % the distance of two parallel lines, which is around the minmum wavelength of stripe 
dst=round(56./(resr));%the distance of two parallel lines for wavelet transformation.
%derived parameters
%Twavmax=Ts1max*1.8; % the maximum wavelength used for wavelet transform
Twavmax = Ts1max; % -> amazing results!
dlmin=Ts1; % the stripes length limit for detecting straight lines.
dlmin2=dlmin/8.; %85m, minimum length of continuous sets of pixels with low gradient
%brightness parameters
brilmt=0.8*max(m0(:));%1800.;
dVqlmt=3;% the maximum difference along stripes 

% Downsizing the image if resolution is better than resr.
nds=round((resr)/(res));
if nds<=1
    flagds=0;  %flag showing downsizing the original image or not, i.e., 1 downsizing, and 0 keeping the original image
    z=m0;
else
    flagds=1; 
[ny,nx]=size(m0);
ixa=[1:nds:nx]';iya=[1:nds:ny]';
[IXAs,IYAs]=meshgrid(ixa,iya);
z=m0(1:nds:end,1:nds:end);
ixa0=[1:nx];iya0=[1:ny];    %store the original index
end
[ny,nx]=size(z);z=double(z);
ixa=[1:nx]';iya=[1:ny]';    
[IXA,IYA]=meshgrid(ixa,iya);

id=find(z>brilmt); idsv1=id;

% % Step 1: Identifying the direction of stripes
% % Criterion: the difference of z along a direction is small, and z > brilmt
% % This step takes 4 seconds of computation time
ds=dsp;
if (opt==1) % direction is known
   angstr=90; %[90,0,-45,45]; % angstr=90; % go back to the version 1, 
else %search for direction
ck1=clock;
% ang=[-89:90]; %the range allowed for angle is (-90, 90]
ang=[];distot=[];
das=[30,10,1];% different step length for line angles in three loops for lesser iterations. %reduce 14sec to 4 sec
for sl=1:length(das) %step loop
    if sl==1;angs=[-89,-60:das(1):60,90];
    else angs=[angsm-das(sl-1)+1:das(sl):angsm+das(sl-1)-1];
    end
    id=find(angs<=-90);angs(id)=angs(id)+180; %same slant angle for a straight line
    id=find(angs>90);angs(id)=angs(id)-180;
    ang=[ang,angs];
    distots=zeros(1,length(angs));
for ia=1:length(angs) % angi=angs(1:end)  %ia=1:1:length(ang) -89:90
    angi=angs(ia);
    [xl,yl,idw,Vq,dsdx]=strline(ixa,iya,z,angi,ds);
    ndlmin=floor(dlmin/dsdx); %number of points within the dlmin
    [lenl,lena]=size(yl);

    %the difference of Vq
    dVq=Vq(2:end)-Vq(1:end-1);
    idd=find(abs(dVq)<=dVqlmt & Vq(1:end-1)>brilmt); %id of difference =0 and Vq < lmt
    iddend=find(mod(idw(idd),lenl)==0);% incase of at the end of a line, delete
    idd(iddend)=[]; % incase of at the end of a line
    %Remove small clusters of data
    Mstr1=zeros(lenl*lena,1);
    Mstr1(idw(idd))=1; % selected data points
    Mstr1 = bwareaopen(Mstr1, ndlmin); % get rid of isolated little clusters of data
    id=find(Mstr1==1);
    %total distances of a constant value;
    distots(ia)=dsdx*length(id);
%     pause
end
distot=[distot,distots];
[maxd,iamd]=max(distot);
iamd=find(distot==maxd);
if(length(iamd)==1)
    angsm=ang(iamd); % the angle with the maximum distance.
else
    display(['Warning: multiple peaks found for loop ',num2str(sl)])
    angsm=mean(ang(iamd));
    das(sl)=floor((max(ang(iamd))-min(ang(iamd)))/2.);%expanding the search range for next loop
end
end
meand=mean(distot);stdd=std(distot);
[maxd,iamd]=max(distot);
maxsign=1; %maximum value is significant; std is not zero %Input control 
if (abs(stdd)<=meand/10. && abs(maxd-meand)<=stdd/10.);maxsign=0;end
if(length(iamd)==1 && maxsign ==1)
    angstr=ang(iamd);
else
    angstr=[90,0,-45,45]; %original grid; consider 4 directions; tends to identify the non-saturation as saturation.
end
ck2=clock;
dt=ck2-ck1;tsec=dt(4)*3600+dt(5)*60+dt(6);
display(['Searching the direction process takes: ',num2str(tsec),' sec'])
end
% % End of Step 1

% Step 2: Preprocessing for wavelet transform.
% Filter out the signal below brilmt, and subtract a constant to make
% the mean to be zero.
id=find(z>brilmt); idsv1=id;
Mc1=z>brilmt; %criterion 1
[ny,nx]=size(z);
mz=mean(z(id));

% Preparing for steps 3 and 4
fc = 5/(2*pi); %fc = centfrq('morl');
angstr=angstr-90; %direction of lines perpendicular to the stripes direction;
id=find(angstr<=-90);angstr(id)=angstr(id)+180; %same slant angle for a straight line
ds=dst;
Mc2ipre=false(ny,nx);nptmax=0;
for ia=1:length(angstr) % angi=angs(1:end)  %ia=1:1:length(ang) -89:90
    ck1=clock;
    angi=angstr(ia);
    
%Step 4: Detecting low brightness gradient areas along the direction of stripes.
    angioth=angi-90; if(angioth<=-90);angioth=angioth+180; end; angi=angioth;
    [xl,yl,idw,Vq,dsdx]=strline(ixa,iya,z,angi,ds);
    [lenl,lena]=size(yl);
    %the difference of Vq
    dVq=Vq(2:end)-Vq(1:end-1);
    idd=find(abs(dVq)<=dVqlmt & Vq(1:end-1)>brilmt); %id of difference =0 and Vq < lmt
    iddend=find(mod(idw(idd),lenl)==0);% incase of at the end of a line, delete
    idd(iddend)=[]; % incase of at the end of a line
    Mstr1=zeros(lenl*lena,1); % to collect selected data points
    Mstr1(idw(idd))=1; % selected data points along the straight lines, 1 means low slope
    nt=round(dlmin2/dsdx); %number of pixels within the dlmin2
    Mstr1 = bwareaopen(Mstr1, nt); % get rid of isolated little clusters of data along lines.
    Mstr1=double(Mstr1);Mstr1=reshape(Mstr1,lenl,lena);
    % interpolation to regular grids on image
    % Note (xp,yp) coordinates are not Cartesian.
    % (xp,yp) coordinate: xp direction along straight lines;
    % yp along y axis when abs(angi)<=45, or along x axis when abs(ang) > 45
    % origin at the origin of the first line with smallest y-intercept.
    if angi==90
        IMc2a= interp2(xl,yl,Mstr1,IXA,IYA,'*linear'); %interpolation along the lines
    elseif (abs(angi)<=45) % use ixa as the input coordinates
        yp0=yl(1,1)-tan(angi*pi/180)*xl(1,1);
        da=yl(1,2)-yl(1,1);
        xp=dsdx*repmat(ixa,1,lena);yp=da*repmat([0:lena-1],lenl,1);
        IXAp=IXA/cos(angi*pi/180);IYAp=IYA-tan(angi*pi/180)*IXA;IYAp=IYAp-yp0;
        IMc2a= interp2(xp',yp',Mstr1',IXAp,IYAp,'*linear');
    else % ang > 45 & <90; use iya as the input coordinates
        yp0=xl(1,1)-yl(1,1)./tan(angi*pi/180);
        dax=xl(1,2)-xl(1,1);
        xp=dsdx*repmat(iya,1,lena);yp=dax*repmat([0:lena-1],lenl,1);
        IXAp=IYA/sin(abs(angi)*pi/180);IYAp=IXA-IYA./tan(angi*pi/180);IYAp=IYAp-yp0;
        IMc2a= interp2(xp',yp',Mstr1',IXAp,IYAp,'*linear');
    end
    id=find(isnan(IMc2a));IMc2a(id)=0;IMc2a=round(IMc2a);
    Mc2a=logical(IMc2a);
    %end of step 4
     
% Step 3: Detecting periodic signals using wavelet transform.
    angi=angstr(ia);
    [xl,yl,idw,Vq,dsdx]=strline(ixa,iya,z,angi,ds);
    ndlmin=floor(dlmin/dsdx); %number of points within the dlmin
    Ts=dsdx;Fs=1./Ts;Trange=[Twavmax 2./Fs];freqrange= 1./Trange;ns=ceil(Trange(1)-Trange(end)+1);
    scalerange = fc./(freqrange*Ts);
    scales = linspace(scalerange(end),scalerange(1),ns); % construct scales at constant pace
    freq = fc./(scales.*Ts);   
    [lenl,lena]=size(yl);
    IMc2str=zeros(lenl,lena);
% Do wavelet transform for each profile
for ka=1:lena % iy=1:ny  %1076
iq=find(idw>=(ka-1)*lenl+1 & idw<=ka*lenl); % find the idw for this particular line 
Vqprof=Vq(iq);
maxVq=max(Vqprof);
if (length(iq)<ndlmin || maxVq < brilmt);continue;end %mininum length of data for detecting stripes
t=(0:length(iq)-1)'*dsdx; % since points within the box are continuous.
prof=zeros(size(Vqprof));Mt=Vqprof>brilmt;
Mt = bwareaopen(Mt, ndlmin); % get rid of isolated little clusters of data
prof(Mt)=Vqprof(Mt)-mz;
%
cwtquadchirp = cwtft({prof,Ts},'wavelet','bump','scales',scales*Ts);
[Sm,id]=max(abs(cwtquadchirp.cfs)); % find index for the dominating wavelength at each location
ix=find(abs(1./cwtquadchirp.frequencies(id))>=Ts1min & abs(1./cwtquadchirp.frequencies(id))<=Ts1max); % ix is the index of id or t, showing location.
idsv2s=idw(iq(ix));
IMc2str(idsv2s)=1;
end %iya
ck2=clock;
dt=ck2-ck1;tsec=dt(4)*3600+dt(5)*60+dt(6);
display(['Wavelet method takes: ',num2str(tsec),' sec',' for angle= ',num2str(angi)])
% interpolation to regular grids on image
if 1 % meshgrid interpolation, fast reduce to 2 seconds from 2 minutes, 
    if angi==90
        IMc2= interp2(xl,yl,IMc2str,IXA,IYA,'*linear'); %interpolation along the lines
    elseif (abs(angi)<=45) % use ixa as the input coordinates
        yp0=yl(1,1)-tan(angi*pi/180)*xl(1,1);
        da=yl(1,2)-yl(1,1);
        xp=dsdx*repmat(ixa,1,lena);yp=da*repmat([0:lena-1],lenl,1);
        IXAp=IXA/cos(angi*pi/180);IYAp=IYA-tan(angi*pi/180)*IXA;IYAp=IYAp-yp0;
        IMc2= interp2(xp',yp',IMc2str',IXAp,IYAp,'*linear');
    else % ang > 45 & <90; use iya as the input coordinates
        yp0=xl(1,1)-yl(1,1)./tan(angi*pi/180);
        dax=xl(1,2)-xl(1,1);
        xp=dsdx*repmat(iya,1,lena);yp=dax*repmat([0:lena-1],lenl,1);
        IXAp=IYA/sin(abs(angi)*pi/180);IYAp=IXA-IYA./tan(angi*pi/180);IYAp=IYAp-yp0;
        IMc2= interp2(xp',yp',IMc2str',IXAp,IYAp,'*linear');
    end
    id=find(isnan(IMc2));IMc2(id)=0;
else  % scatter points interpolation, slow
    F= scatteredInterpolant(xl(:),yl(:),IMc2str(:)); %interpolation, accept input grid as scatter points
    IMc2 =F(IXA,IYA);
end
IMc2=round(IMc2);
%End of step 3

Mc2i=logical(IMc2)&Mc2a; % combination os steps 3 and 4

%Choose only one direction that gives more points.
npt=sum(Mc2i(:));
if npt>nptmax; 
    Mc2=Mc2i;Mc2ipre=Mc2;nptmax=npt;
else
    Mc2=Mc2ipre;
end
end %ia=1:length(angstr)

M=Mc1&Mc2;

if flagds==0 % No need of downsizing 
    npt=sum(sum((z~=0)));nptsat=sum(sum(M));
    if npt==0; ratio=0;else ratio=nptsat/npt*100;end
else % size back to original image, compatible with datamask, 
    ck4=clock;
    mssd=double(M);
    mssdl= interp2(IXAs,IYAs,mssd,ixa0,iya0','*linear',0); %linear produce NaNs for domain outside.  A(isnan(A)) = 0 ;
    M=logical(round(mssdl));
    ck5=clock;
    dt=ck5-ck4;tsec=dt(4)*3600+dt(5)*60+dt(6);
    display(['Size recovering process takes: ',num2str(tsec),' sec'])
    npt=sum(sum((m0~=0)));nptsat=sum(sum(M));
    if npt==0; ratio=0;else ratio=nptsat/npt*100;end
end

cksatu2=clock;
dt=cksatu2-cksatu1;tsec=dt(4)*3600+dt(5)*60+dt(6);
display(['Saturation filtering process takes: ',num2str(tsec),' sec'])

end %end of function

function [xl,yl,idw,Vq,dsdx]=strline(ixa,iya,z,angi,ds)
    % get coordinates of straight lines with slant angle as angi 
    % angi, the slant angle of straight lines, in degrees (-89, 90];
    % ds, the distance of two parallel lines in pixels;
    % ixa, the coordinates of the image area along x-axis;
    % iya, the coordinates of the image area along y-axis;
    % z, the brightness of the image, 2-D matrix
    % Output parameters:
    % xl, x-coordinates of the straight lines; 
    % yl, y-coordinates of the straight lines;
    %     note: xl,yl has the same coordinate system as the input, ixa, iya.
    % idw, index of the xl, yl that are within the image area;
    % Vq, the interpolated brightness on the straight lines;
    % dsdx, the distance between two sequential points on a straight line
    [ny,nx]=size(z);
    if angi==90
        lena=length(ixa(1:ds:end)');        
        xl= repmat(ixa(1:ds:end)',ny,1); 
        yl=repmat(iya,1,lena);
        idw=find(xl>=ixa(1) & xl<=ixa(end)); 
        Xq=xl(idw);Yq=yl(idw);Vq=z(:,1:ds:end);Vq=Vq(:);
        dsdx=1;
    else
    if angi>=0 ; xc=ixa(1);else xc=ixa(end); end
    a=[-tan(angi*pi/180)*ixa(1:nx-1:end)',iya(1:ny-1:end)'-tan(angi*pi/180)*xc];% points (xc,iya) (ixa,0)
    da=ds/cos(angi*pi/180);
    a=min(a):da:max(a);
    lena=length(a); 
    if (abs(angi)<=45) % use ixa as the input coordinates
        bx=tan(angi*pi/180)*ixa;
        [A,BX]=meshgrid(a,bx);
        yl=A+BX; %xl,yl the coordinates of the straight line
        xl= repmat(ixa,1,lena); 
        [idw]=find(yl>=iya(1) & yl<=iya(end)); % array inside of the box
        Xq=xl(idw);Yq=yl(idw);
    else % ang > 45; use iya as the input coordinates
        yl=repmat(iya,1,lena);A=repmat(a,ny,1);b=tan(angi*pi/180);
        xl=(yl-A)/b;
        idw=find(xl>=ixa(1) & xl<=ixa(end)); 
        Xq=xl(idw);Yq=yl(idw);
    end
    dsdx=abs(xl(2)-xl(1))/cos(angi*pi/180); % distance of two sequential points on the line
    [IXA,IYA]=meshgrid(ixa,iya);
    Vq = interp2(IXA,IYA,z,Xq,Yq,'*linear'); %interpolation along the lines
    end % ang == 90
    
end


