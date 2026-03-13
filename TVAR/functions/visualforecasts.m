clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% strFolder = 'Y:\Users\146431\pier\december2012\TARCOMPARISON\';
% strFolder = 'D:\My_Documents\Density_forecasts\Code\TVAR\files\';
strFolder = 'I:\FCI';

% addpath(strcat(strFolder,'functions'));

% Directories with data and forecasts (one for each model)
dirname     = strcat(strFolder,'data',filesep,'databenchmark',filesep); 
fdirnametar = strcat(strFolder,'forecasts',filesep,'TAR',filesep);
fdirnamevar = strcat(strFolder,'forecasts',filesep,'benchmark',filesep);
fdirnamebasic = strcat(strFolder,'forecasts',filesep,'basic',filesep); 
fdirnametvtp = strcat(strFolder,'tvtp',filesep,'forecasts',filesep); 


% Number of the first/last data file and number of retained draws:
minfile = 0;
maxfile = 354;
S       = 5000;

% Number of observations used in initial estimation, horizons, density
% percentiles to be saved: 
T0       = 121 ;
horizons = [1 6 12];
percentiles = [10 25 50 75 90];
% NB must be the SAME FOR ALL MODELS

% Load data:
tempname = strcat(dirname,'data',num2str(maxfile),'.mat');
load(tempname);
datavar = dataout;
datatar = dataout;
databasic = dataout(:,1:3);
datatvtp = dataout(:,1:3);


% NB the data can DIFFER ACROSS MODELS. Eg datatar could be set by:
% tempname = (...);
% load(tempname);
% datatar = dataout;

% Define labels:
models    = {'VAR ', 'TAR', 'BASIC'};
vnamesvar = {'ip ', 'r ', 'cpi ', 'fci '};
vnamestar = vnamesvar ;
vnamesbasic = {'ip ', 'r ', 'cpi '};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialise outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = rows(datavar);
N = cols(datavar);
T2 = maxfile + 1 - max(horizons);
% T2 is the # of usable forecasts. If minfile=0 we need to add +1, because
% we have [forecast0 .. forecastK] = (K+1) forecasts.


dates=1983+(2/12):1/12:2009+(1/12);
[tesmp]=xlsread('dates.xls');
dates=datestr(tesmp,2);

set(gca,'nextplot','replacechildren');

%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mm=1;
for file = 287 : 311  %2007-01-01 to 2009-01-01


   
    % Load draws
    %----------------------------------------------------------------------
	fnameTAR=strcat(fdirnametar,'forecast',num2str(file),'.mat');
    fnameVAR=strcat(fdirnamevar,'forecast',num2str(file),'.mat');
    fnamebasic=strcat(fdirnamebasic,'forecast',num2str(file),'.mat');
        fnametvtp=strcat(fdirnametvtp,'forecast',num2str(file),'.mat');



    load(fnameTAR);
    ftar = fsave;

    load(fnameVAR);
    fvar = fsave;
    
    load(fnamebasic);
    fbasic = fsave;
    
    load(fnametvtp);
    ftvtp = fsave;
    
	% Split data into estimation sample and observations for fc/density evaluation
    %----------------------------------------------------------------------
    samplevar = datavar(1:T0+file, :);
    sampletar = datatar(1:T0+file, :);
    samplebasic = databasic(1:T0+file, :);
        sampletvtp = datatvtp(1:T0+file, :);


    
    obsvar = datavar(T0+file+1 : T0+file+max(horizons), :) ;
    obstar = datatar(T0+file+1 : T0+file+max(horizons), :) ; 
    obsbasic = databasic(T0+file+1 : T0+file+max(horizons), :) ; 
        obstvtp = datatvtp(T0+file+1 : T0+file+max(horizons), :) ; 

    
    
    % Evaluate point forecasts
    %----------------------------------------------------------------------
    pfvar   = squeeze(mean(fvar,1));
    pftar   = squeeze(mean(ftar,1));
    pfbasic   = squeeze(mean(fbasic,1));

% figure(file)    

% fig=figure('Name',char(dates(file)),'NumberTitle','off');
fig=figure;
v=1; %IP
subplot(3,4,1)
plotfan(fvar(:,:,v)')
hold on
plot(obsvar(:,v),'b');
axis tight
title('BVAR');
ylabel('IP')
temp=prctile(fvar(:,:,v),[ 5 95],1);
text(-0.01,max(max(temp))+1,(dates(file,:)),'FontSize',14)

subplot(3,4,3)
plotfan(ftar(:,:,v)')
hold on
plot(obstar(:,v),'b');
axis tight

title('TAR');


subplot(3,4,2)
plotfan(fbasic(:,:,v)')
hold on
plot(obsbasic(:,v),'b');
axis tight

title('BVAR BASIC');

subplot(3,4,4)
plotfan(ftvtp(:,:,v)')
hold on
plot(obstvtp(:,v),'b');
axis tight
title('MSVAR TVTP');





v=2; %R
subplot(3,4,5)
plotfan(fvar(:,:,v)')
hold on
plot(obsvar(:,v),'b');
axis tight
title('BVAR');
ylabel('R')
text(min(min(fvar(:,:,1)))-0.01,max(max(fvar(:,:,1)))+0.01,char(dates(file)))

subplot(3,4,7)
plotfan(ftar(:,:,v)')
hold on
plot(obstar(:,v),'b');
axis tight

title('TAR');


subplot(3,4,6)
plotfan(fbasic(:,:,v)')
hold on
plot(obsbasic(:,v),'b');
axis tight

title('BVAR BASIC');

subplot(3,4,8)
plotfan(ftvtp(:,:,v)')
hold on
plot(obstvtp(:,v),'b');
axis tight
title('MSVAR TVTP');









v=3; %PI
subplot(3,4,9)
plotfan(fvar(:,:,v)')
hold on
plot(obsvar(:,v),'b');
axis tight
title('BVAR');
ylabel('\pi')
text(min(min(fvar(:,:,1)))-0.01,max(max(fvar(:,:,1)))+0.01,char(dates(file)))

subplot(3,4,11)
plotfan(ftar(:,:,v)')
hold on
plot(obstar(:,v),'b');
axis tight

title('TAR');


subplot(3,4,10)
plotfan(fbasic(:,:,v)')
hold on
plot(obsbasic(:,v),'b');
axis tight

title('BVAR BASIC');

subplot(3,4,12)
plotfan(ftvtp(:,:,v)')
hold on
plot(obstvtp(:,v),'b');
axis tight
title('MSVAR TVTP');


















 fmat(mm)=getframe(fig);
mm=mm+1;    
    disp(file);
end

movie2avi(fmat, 'forecasts.avi', 'compression', 'None','fps',2);

close all
[h, w, p] = size(fmat(1).cdata); % use 1st frame to get dimensions
hf = figure; 
% resize figure based on frame's w x h, and place at (150, 150)
set(hf, 'position', [150 150 w h]);
axis off
movie(hf,fmat);



