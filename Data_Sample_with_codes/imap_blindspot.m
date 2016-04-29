% This demo script reproduces the analysis reported in the user guidebook
% Junpeng Lao, 2015 March, University of Fribourg
%--------------------------------------------------------------------------
% Copyright (C) iMap Team 2015
%% Load data
clear all
clc
cd ('Data_sample_blindspot')
filename  = dir('data*.mat');
% create condition table
sbj       = repmat([1:30],1,4)';
group     = repmat([ones(1,15) ones(1,15)*2],1,4)';
blinkspot = repmat(1:4,30,1);
blinkspot = blinkspot(:);
Tbl       = dataset(sbj,group,blinkspot);
% deg0cau = 1:15;
% deg0as  = 16:30;
% deg2cau = 31:45;
% deg2as  = 46:60;
% deg5cau = 61:75;
% deg5as  = 76:90;
% deg8cau = 91:105;
% deg8as  = 106:120;

% parameters for smoothing
ySize        = 382;
xSize        = 390;
smoothingpic = 10;
[x, y]       = meshgrid(-floor(xSize/2)+.5:floor(xSize/2)-.5, -floor(ySize/2)+.5:floor(ySize/2)-.5);
gaussienne   = exp(- (x .^2 / smoothingpic ^2) - (y .^2 / smoothingpic ^2));
gaussienne   = (gaussienne - min(gaussienne(:))) / (max(gaussienne(:)) - min(gaussienne(:)));
f_fil        = fft2(gaussienne);

Nitem        = length(filename);
rawmapMat    = zeros(Nitem, ySize, xSize);
fixmapMat    = zeros(Nitem, ySize, xSize);
stDur        = zeros(Nitem,1);
Ntall        = zeros(Nitem,1);
for item = 1:Nitem
    load(['data' num2str(item) '.mat'])
    Nfix           = size(summary, 1);
    Trials         = unique(summary(:, 4));
    Ntall(item)    = length(Trials);
    % condition based
    coordX         = round(summary(:, 2));%switch matrix coordinate here
    coordY         = round(summary(:, 1));
    intv           = summary(:, 3);
    indx1          = coordX>0 & coordY>0 & coordX<ySize & coordY<xSize;
    rawmap         = full(sparse(coordX(indx1),coordY(indx1),intv(indx1),ySize,xSize));
    f_mat          = fft2(rawmap); % 2D fourrier transform on the points matrix
    filtered_mat   = f_mat .* f_fil;
    smoothpic           = real(fftshift(ifft2(filtered_mat)));
    fixmapMat(item,:,:) = (smoothpic-mean(smoothpic(:))) ./ std(smoothpic(:));%smoothpic;% ./sum(durind); % 
    rawmapMat(item,:,:) = rawmap;
    stDur(item)         = sum(indx1);
end
Tbl.stDur                 = stDur;
Tbl.blinkspot             = nominal(Tbl.blinkspot);
Tbl.group                 = nominal(Tbl.group);
Tbl.group(Tbl.group=='1') = 'WC';
Tbl.group(Tbl.group=='2') = 'EA';
Tbl.group                 = nominal(Tbl.group);

%% rescale
scale              = 150 / mean([xSize,ySize]);
[ySize2, xSize2]   = size(imresize(ones(ySize,xSize), scale, 'nearest'));
fixmapMat2         = zeros(size(Tbl, 1), ySize2, xSize2);
for it = 1:size(Tbl,1)
    fixmapMat2(it,:,:) = imresize(squeeze(fixmapMat(it,:,:)), scale, 'nearest');
end

%% Mean fixation intensity map by condition
figure('NumberTitle','off','Name','Mean fixation bias');
race       = unique(Tbl.group);
condi      = unique(Tbl.blinkspot);
for ig = 1:length(race)
    for ipp = 1:length(condi)
        indxtmp = (Tbl.group == race(ig)) & (Tbl.blinkspot == condi(ipp));
        TempMap = squeeze(mean(fixmapMat2(indxtmp,:,:), 1));
        subplot(2,4,(ig-1)*4+ipp)
        imagesc(TempMap);colorbar
        axis square, axis off,
        title(char(race(ig)))
    end
end
%
figure('NumberTitle','off','Name','all fixation bias');
subplot(1,2,1)
imagesc(squeeze(mean(fixmapMat2,1)));
axis('equal','off');
subplot(1,2,2)
masktmp=squeeze(mean(fixmapMat2,1))>.01;
imagesc(masktmp);
axis('equal','off');
%% imapLMM
% tic
% opt.singlepredi=1;
% [LMMmap,lmexample]=imapLMM(fixmapMat2,Tbl,masktmp,opt, ...
%    'PixelIntensity ~ group + blinkspot +  group:blinkspot + (1|sbj)','DummyVarCoding','effect');
% toc
% %% imapLMM_result
% %% plot model fitting
% close all
% opt1.type='model';
% % perform contrast
% [StatMap]=imapLMMcontrast(LMMmap,opt1);
% % output figure;
% imapLMMdisplay(StatMap,0)
%% Single trial
Trialall    = sum(Ntall);
rawmapMatST = zeros(Trialall, ySize, xSize);
fixmapMatST = zeros(Trialall, ySize, xSize);

groupST     = zeros(Trialall,1);
sbjST       = zeros(Trialall,1);
blinkspotST = zeros(Trialall,1);
stDurST     = zeros(Trialall,1);
ii = 0;
for item = 1:Nitem
    load(['data' num2str(item) '.mat'])
    Nfix    = size(summary,1);
    Trials  = summary(:,4);
    UTrials = unique(Trials);
    Nt      = length(UTrials);
    % condition based
    fixmaptmp = zeros(Nt, ySize, xSize);
    rawmaptmp = zeros(Nt, ySize, xSize);
    for it = 1:Nt
        ii     = ii+1;
        coordX = round(summary(Trials == UTrials(it), 2));%switch matrix coordinate here
        coordY = round(summary(Trials == UTrials(it), 1));
        intv   =       summary(Trials == UTrials(it), 3);
        
        indx1  = coordX>0 & coordY>0 & coordX<ySize & coordY<xSize;
        rawmap = full(sparse(coordX(indx1),coordY(indx1),intv(indx1),ySize,xSize));
        f_mat  = fft2(rawmap); % 2D fourrier transform on the points matrix
        filtered_mat        = f_mat .* f_fil;
        smoothpic           = real(fftshift(ifft2(filtered_mat)));
        fixmaptmp(it,:,:)   = (smoothpic-mean(smoothpic(:)))./std(smoothpic(:));%smoothpic./sum(indx1); % 
        fixmapMatST(ii,:,:) = fixmaptmp(it,:,:);
        rawmaptmp(it,:,:)   = rawmap;
        rawmapMatST(ii,:,:) = rawmaptmp(it,:,:);
        
        groupST(ii)         = Tbl.group(item);
        sbjST(ii)           = Tbl.sbj(item);
        blinkspotST(ii)     = Tbl.blinkspot(item);
        stDurST(ii)         = sum(indx1);
    end
    if Nt<10
        fixmapMat(item,:,:) = mean(fixmaptmp,1);
        rawmapMat(item,:,:) = mean(rawmaptmp,1);
    else
        fixmapMat(item,:,:) = trimmean(fixmaptmp,30);
        rawmapMat(item,:,:) = trimmean(rawmaptmp,30);
    end
end
TblST                             = dataset(sbjST,groupST,blinkspotST,stDurST);
TblST.blinkspotST                 = nominal(TblST.blinkspotST);
TblST.groupST                     = nominal(TblST.groupST);
TblST.groupST(TblST.groupST=='1') = 'WC';
TblST.groupST(TblST.groupST=='2') = 'EA';
TblST.groupST                     = nominal(TblST.groupST);

%% rescale
% scale=150/mean([xSize,ySize]);
% [ySize2,xSize2]=size(imresize(ones(ySize,xSize),scale,'nearest'));
% fixmapMat2ST=zeros(size(TblST,1),ySize2,xSize2);
% for it=1:size(TblST,1)
%     fixmapMat2ST(it,:,:)=imresize(squeeze(fixmapMatST(it,:,:)),scale,'nearest');
% end
%% mean map
figure('NumberTitle','off','Name','Mean fixation bias');
race        = unique(TblST.groupST);
condi       = unique(TblST.blinkspotST);
ii = 0;
for ig = 1:length(race)
    for ipp = 1:length(condi)
        indxtmp = (TblST.groupST == race(ig)) & (TblST.blinkspotST == condi(ipp));
        TempMap = squeeze(mean(fixmapMatST(indxtmp, :, :), 1));
        % subplot(2,4,(ig-1)*4+ipp)
        subplot(4,2,ipp*2-(2-ig))
        ii              = ii+1;
        betanew(ii,:,:) = TempMap;
        imagesc(TempMap);colorbar
        axis square, axis off,
        title(char(race(ig)))
    end
end
%
figure('NumberTitle','off','Name','all fixation bias');
subplot(1,2,1)
imagesc(squeeze(mean(fixmapMatST,1)));
axis('equal','off');
subplot(1,2,2)
masktmpST=squeeze(mean(fixmapMatST,1))>.0045;
imagesc(masktmpST);
axis('equal','off');

%% imapLMM
tic
opt.singlepredi  = 1;
% opt.parallelname = 'grid2';
% [LMMmap,lmexample] = imapLMM(fixmapMat2,Tbl,masktmp,opt, ...
%    'PixelIntensity ~ group + blinkspot +  group:blinkspot + (stDur|sbj)','DummyVarCoding','effect');
[LMMmap, lmexample] = imapLMM(fixmapMatST,TblST,masktmpST,opt, ...
   'PixelIntensity ~ groupST + blinkspotST +  groupST:blinkspotST + (1|sbjST)',...
   'DummyVarCoding','effect');
toc
%% plot fixed effec(anova result)
% close all
opt       = struct;% clear structure
opt.type  = 'fixed';
opt.alpha = 0.05;

mccopt             = struct;
mccopt.methods     = 'bootstrap';
mccopt.bootgroup   = {'groupST'};
% mccopt.methods     = 'FDR';
mccopt.nboot       = 1000;
% mccopt.permute     = 1;
mccopt.bootopt     = 1;
% mccopt.sigma       = smoothingpic*scale;
mccopt.tfce        = 0;

% perform contrast
[StatMap]   = imapLMMcontrast(LMMmap, opt);
[StatMap_c] = imapLMMmcc(StatMap, LMMmap, mccopt, fixmapMatST);
% output figure;
imapLMMdisplay(StatMap_c,1)
%% Compute linear contrast (reproduce figure 2 as in the orignial paper)
% close all
opt       = struct;% clear structure
opt.type  = 'predictor beta';
opt.alpha = 0.05;
opt.c     = {[-1 0 0 0 1 0 0 0]; ...
             [0 -1 0 0 0 1 0 0]; ...
             [0 0 -1 0 0 0 1 0]; ...
             [0 0 0 -1 0 0 0 1]; ...
             [1 0 0 -1 0 0 0 0]; ...
             [0 0 0 0 1 0 0 -1]};
opt.name  = {'WC-EA NV'; 'WC-EA 2dg'; 'WC-EA 5dg'; 'WC-EA 8dg'; 'WC NV-8dg'; 'EA NV-8dg'};
% opt.c       = limo_OrthogContrasts([3,2]);
% opt.name    = {'Spotlight';'Position';'Interaction'};
% opt.h       = {[0.005],[0.005],[0.005],[0.005],[0.005],[0.005]};
% opt.onetail = '>';

mccopt             = struct;
% mccopt.methods     = 'Randomfield';
mccopt.methods     = 'bootstrap';
mccopt.bootgroup   = {'groupST'};
mccopt.nboot       = 1000;
% mccopt.permute     = 1;
mccopt.bootopt     = 1;
% mccopt.sigma       = smoothingpic*scale;
mccopt.tfce        = 0;

% perform contrast
[StatMap]   = imapLMMcontrast(LMMmap,opt);
[StatMap_c] = imapLMMmcc(StatMap,LMMmap,mccopt,fixmapMatST);

% output figure;
imapLMMdisplay(StatMap_c,1)
%% Single predictor (above chance fixate)
opt       = struct;% clear structure
opt.type  = 'predictor beta';
opt.alpha = 0.05;
opt.c     = {[1 0 0 0 0 0 0 0]; ...
             [0 1 0 0 0 0 0 0]; ...
             [0 0 1 0 0 0 0 0]; ...
             [0 0 0 1 0 0 0 0]; ...
             [0 0 0 0 1 0 0 0]; ...
             [0 0 0 0 0 1 0 0]; ...
             [0 0 0 0 0 0 1 0]; ...
             [0 0 0 0 0 0 0 1]};
opt.name  = {'EA-NV'; 'EA-2deg'; 'EA-5deg'; 'EA-8deg'; 'WC-NV'; 'WC-2deg'; 'WC-5deg'; 'WC-8deg'};
h0        = mean(fixmapMatST(repmat(masktmpST, [size(fixmapMatST,1),1,1]) == 1));
opt.h     = {h0, h0, h0, h0, h0, h0, h0, h0};
% opt.c     = limo_OrthogContrasts([3,2]);
% opt.name  = {'Spotlight';'Position';'Interaction'};
% opt.h     = {[0.005],[0.005],[0.005],[0.005],[0.005],[0.005]};
opt.onetail = '>';

mccopt             = struct;
% mccopt.methods     = 'Randomfield';
mccopt.methods     = 'bootstrap';
mccopt.bootgroup   = {'groupST'};
mccopt.nboot       = 1000;
% mccopt.permute     = 1;
mccopt.bootopt     = 1;
% mccopt.sigma       = smoothingpic*scale;
mccopt.tfce        = 0;

% perform contrast
[StatMap]   = imapLMMcontrast(LMMmap,opt);
[StatMap_c] = imapLMMmcc(StatMap,LMMmap,mccopt,fixmapMatST);

% output figure;
imapLMMdisplay(StatMap_c,1,[],'parula')
%% bootCI for fixed effect beta
% size(StatMap_c.betaCI)
% size(StatMap_c.resampMat.resBeta)
cmax   = max(StatMap_c.betaCI(:));
cmin   = min(StatMap_c.betaCI(:));
alpha  = 0.05;
[nboot,Nc,xSize,ySize] = size(StatMap_c.resampMat.resBeta);

for ic=1:Nc
    betaORI = squeeze(StatMap_c.beta(ic,:,:));
    figure('Numbertitle','off','Name',StatMap_c.label{ic});
    betaCI1 = squeeze(StatMap_c.betaCI(ic,1,:,:));
    betaCI2 = squeeze(StatMap_c.betaCI(ic,2,:,:));
    subplot(2,3,1)
    imagesc(betaCI1,[cmin cmax]);axis('equal','off')
    title('lower CI')
    subplot(2,3,4)
    imagesc(betaCI2,[cmin cmax]);axis('equal','off')
    title('upper CI')
    % compute Bootstrap CI
    betatmp  = squeeze(StatMap_c.resampMat.resBeta(:,ic,:,:));
    betasort = sort(betatmp(:,:),1);
    betaCI1b = reshape(betasort(nboot*alpha,:),     [xSize,ySize]) + betaORI;
    betaCI2b = reshape(betasort(nboot*(1-alpha),:), [xSize,ySize]) + betaORI;
    
    subplot(2,3,2)
    imagesc(betaCI1b,[cmin cmax]);axis('equal','off')
    title('lower CI (bs)')
    subplot(2,3,5)
    imagesc(betaCI2b,[cmin cmax]);axis('equal','off')
    title('upper CI (bs)')

    subplot(2,3,3)
    imagesc(betaCI1b - betaCI1);colorbar;axis('equal','off')
    title('diff CI (bs-pm)')
    subplot(2,3,6)
    imagesc(betaCI2b - betaCI2);colorbar;axis('equal','off')
    title('diff CI (bs-pm)')
end
