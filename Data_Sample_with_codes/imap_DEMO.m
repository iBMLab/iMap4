% This demo script reproduces the analysis reported in the user guidebook
% Junpeng Lao, 2015 March, University of Fribourg
%--------------------------------------------------------------------------
% Copyright (C) iMap Team 2015
%% load data file
p = userpath;
addpath(genpath([ p(1:end-1) '/Apps/iMAP']));
clear all
cd('./Data_sample_DEMO/')
datafile = 'data_subset.txt';
% data structure
% column 1 :  subject
% column 2 :  TRIAL
% column 3 :  xfix
% column 4 :  yfix
% column 5 :  FIX_START
% column 6 :  FIX_END
% column 7 :  fixdur
% column 8 :  Spotlight
% column 9 :  Rating
% column 10:  IMAGE
% column 11:  Position
% column 12:  EYE_USED

filetxt     = fopen(datafile);
data        = textscan(filetxt,'%s%n%n%n%n%n%n%s%n%s%s%s%s');
subject     = (data{1});
TRIAL       = (data{2});
xfix        = (data{3});
yfix        = (data{4});
FIX_START   = (data{5});
FIX_END     = (data{6});
fixdur      = (data{7});
Spotlight   = (data{8});
Rating      = (data{9});
stimli      = (data{10});
Position    = (data{11});
%% experiment parameters
dist        = 83;
yScreen     = 30;
xScreen     = 37;
ySize       = 1024;
xSize       = 1280;
% rescale
scale       = 1/5;
[ySize2, xSize2] = size(imresize(ones(ySize,xSize),scale,'nearest'));
% create Gaussian for smoothing
smoothingpic     = 20;
[x, y]           = meshgrid(-floor(xSize/2)+.5:floor(xSize/2)-.5,...
                            -floor(ySize/2)+.5:floor(ySize/2)-.5);
gaussienne       = exp(- (x .^2 / smoothingpic ^2) - (y .^2 / smoothingpic ^2));
gaussienne       = (gaussienne - min(gaussienne(:))) ...
                 / (max(gaussienne(:)) - min(gaussienne(:)));
f_fil            = fft2(gaussienne);

subjlist         = unique(subject);
Ns               = length(subjlist);% number of subject
Nc               = 3;%number of condition (Natural viewing, 100, 300)
Nsti             = length(unique(stimli));% number of image
%% prepare fixation matrix
SLname           = unique(Spotlight);
Nc1              = length(SLname);
PSname           = unique(Position);
Nc2              = length(PSname);
Nitem            = Nc1*Nc2*Ns;
FixMap           = zeros(Nitem, ySize2, xSize2);
RawMap           = zeros(Nitem, ySize2, xSize2);
indx             = zeros(Nitem, 1);
ii               = 0;
for is = 1:Ns
    disp(is)
    % read individual data
    indSbj       = strcmp(subject, subjlist{is});
    seyex        = xfix(indSbj);
    seyey        = yfix(indSbj);
    sfixs        = FIX_START(indSbj);
    sfixd        = fixdur(indSbj);
    strial       = TRIAL(indSbj);
    scondi       = Spotlight(indSbj);
    srate        = Rating(indSbj);
    simag        = stimli(indSbj);
    spost        = Position(indSbj);
    tbl          = dataset(seyex,seyey,sfixs,sfixd,strial,scondi,srate,simag,spost);
    % figure;plot(seyex,ySize-seyey,'.r');axis([0 xSize 0 ySize])
    %% convolution with Gaussian kernel
    for ic1 = 1:Nc1
        for ic2 = 1:Nc2
            %%
            ii       = ii+1;
            selected = strcmp(scondi, SLname(ic1)) & strcmp(spost, PSname(ic2));
            if sum(selected) ~= 0
                trialid       = unique(strial(selected));
                isfixmap      = zeros(size(trialid,1), ySize2, xSize2);
                israwmap      = isfixmap;
                stDur         = zeros(size(trialid,1), 1);
                stRate        = stDur;
                for it = 1:size(trialid,1)
                    % fixation matrix
                    selected2 = find(strial == trialid(it));
                    seyext    = seyex(selected2);
                    seyeyt    = seyey(selected2);
                    intv      = sfixd(selected2);
                    
                    rawmap    = zeros(ySize, xSize);
                    coordX    = round(seyeyt);%switch matrix coordinate here
                    coordY    = round(seyext);
                    indx1=coordX>0 & coordY>0 & coordX<ySize & coordY<xSize;
                    if isempty(indx1)
                        indx(ii) = 1;
                    else
                        rawmap       = full(sparse(coordX(indx1), coordY(indx1), intv(indx1), ySize, xSize));
                        f_mat        = fft2(rawmap); % 2D fourrier transform on the points matrix
                        filtered_mat = f_mat .* f_fil;
                        
                        smoothpic        = real(fftshift(ifft2(filtered_mat)));
                        isfixmap(it,:,:) = imresize(smoothpic, scale, 'nearest');
                        israwmap(it,:,:) = imresize(rawmap,    scale, 'nearest');
                        stDur (it)       = sum(intv(indx1));
                        stRate(it)       = srate(selected2(1));
                    end
                end
                Subject{ii,1}            = subjlist{is};
                spotlight(ii,1)          = SLname(ic1);
                position(ii,1)           = PSname(ic2);
                rating(ii,1)             = nanmean(stRate);
                fixDur(ii,1)             = nanmean(stDur);
                FixMap(ii,:,:)           = nanmean(isfixmap,1);
                RawMap(ii,:,:)           = nanmean(israwmap,1);
            else
                indx(ii) = 1;
            end
        end
    end
end
%% matrix for inputing into imapLMM
PredictorM             = dataset(Subject,spotlight,position,rating,fixDur);
indx2                  = isnan(FixMap(:,1,1));
FixMap    (logical(indx+indx2),:,:)   = [];
RawMap    (logical(indx+indx2),:,:)   = [];
PredictorM(logical(indx+indx2),:)     = [];
Mask                 = squeeze(nanmean(FixMap))>(min(fixdur)/2);
PredictorM.Subject   = nominal(PredictorM.Subject);
PredictorM.spotlight = nominal(PredictorM.spotlight);
PredictorM.position  = nominal(PredictorM.position);
%% mean fixation map
figure('NumberTitle','off','Name','Mean fixation bias');
Spotlightsize  = unique(PredictorM.spotlight);
Position       = unique(PredictorM.position);
ord            = [1 3 5 2 4 6];
io             = 0;
for isp = 1:3
    for ipp = 1:2
        indxtmp    = PredictorM.spotlight == Spotlightsize(isp) ...
                   & PredictorM.position  == Position(ipp);
        TempMap    = squeeze(nanmean(FixMap(indxtmp,:,:), 1));
        
        subplot(2,3,(ipp-1)*3+isp);
        imagesc(TempMap);
        axis off,axis equal;
        title(char([Position(ipp) Spotlightsize(isp)]));
        io = io+1;
    end
end
% figure;imagesc(Mask)
% title('mask')
% axis off,axis equal
%% Linear Mixed Modeling with imapLMM
tic
opt.singlepredi    = 1;
[LMMmap,lmexample] = imapLMM(FixMap,PredictorM,Mask,opt, ...
    'PixelIntensity ~ spotlight + position + spotlight:position + (fixDur|Subject)', ...
    'DummyVarCoding','effect','FitMethod','REML');
toc
%% LMMmap structure
% imapLMM will output a structure called LMMmap with all the information
% from the linear mixed model fitting:
% LMMmap =
%
%                    runopt: [1x1 struct]
%              VariableInfo: [6x4 dataset]
%                 Variables: [118x6 dataset]
%                 FitMethod: 'REML'
%                   Formula: [1x1 classreg.regr.LinearMixedFormula]
%                    modelX: [118x6 double]
%                FitOptions: {'DummyVarCoding'  'effect' 'Fitmethod' 'REML'}
%                  modelDFE: 112
%          CoefficientNames: {1x6 cell}
%                     Anova: [1x1 struct]
%                SinglePred: [1x1 struct]
%             RandomEffects: [1x1 struct]
%     CoefficientCovariance: [4-D double]
%                       MSE: [256x320 double]
%                       SSE: [256x320 double]
%                       SST: [256x320 double]
%                       SSR: [256x320 double]
%                  Rsquared: [2x256x320 double]
%            ModelCriterion: [4x256x320 double]
%              Coefficients: [4-D double]
%% show result
%% plot model fitting
close all
opt1.type  = 'model';
% perform contrast
[StatMapm] = imapLMMcontrast(LMMmap,opt1);
% output figure;
imapLMMdisplay(StatMapm,0,'front.jpg');
%% plot fixed effec(anova result)
% close all
opt         = struct;% clear structure
opt.type    = 'fixed';
opt.alpha   = 0.05;
% perform contrast
[StatMap]   = imapLMMcontrast(LMMmap,opt);

mccopt          = struct;
mccopt.methods  = 'bootstrap';
mccopt.nboot    = 1000;
mccopt.bootopt  = 1;
mccopt.tfce     = 0;
% perform multiple comparison correction
[StatMap_c]     = imapLMMmcc(StatMap,LMMmap,mccopt,FixMap);
% output figure;
imapLMMdisplay(StatMap_c,1,'front.jpg',[])
%% post-hoc
[PostHoc]       = imapLMMposthoc(StatMap_c,RawMap,LMMmap,'mean');
%% Linear Contrast
% close all
opt       = struct;% clear structure
opt.type  = 'predictor beta';
opt.alpha = 2e-07;
opt.c{1}  = [1 0 0 0 -1 0];
opt.name  = {'back 100-NV'};
% perform contrast
[StatMap] = imapLMMcontrast(LMMmap,opt);

mccopt          = struct;
mccopt.methods  = 'bootstrap';
mccopt.nboot    = 1000;
mccopt.bootopt  = 1;
mccopt.tfce     = 0;
% multiple comparison correction
[StatMap_c]     = imapLMMmcc(StatMap,LMMmap,mccopt,FixMap);
% output figure;
imapLMMdisplay(StatMap_c,0,'back.jpg')
%% One-tail t-test
% close all
opt          = struct;% clear structure
opt.type     = 'predictor beta';
opt.alpha    = 0.05;
opt.onetail  = '>';
opt.h        = mean(mean(FixMap(:,Mask)));
opt.c{1}     = [0 1 0 0 0 0];
opt.name     = {'front 100'};
% perform contrast
[StatMap]    = imapLMMcontrast(LMMmap,opt);

mccopt          = struct;
mccopt.methods  = 'bootstrap';
mccopt.nboot    = 1000;
mccopt.bootopt  = 1;
mccopt.tfce     = 0;
% multiple comparison correction
[StatMap_c]     = imapLMMmcc(StatMap,LMMmap,mccopt,FixMap);
% output figure;
imapLMMdisplay(StatMap_c,0,'front.jpg')
