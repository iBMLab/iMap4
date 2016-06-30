% the following script is a demonstration of the Single-Trial estimation
% method for iMap4.
% In this hypothetical experiment, two groups of subjects participated in a
% free viewing experiment of face recognition. We introduced a main effect
% between groups that control subject display a triangle pattern whereas
% the patient group only look at the month. However, there is an effect of
% the eye region for the patient group only: if they fixated on the eye
% they gave more accurate response (Interaction between group and
% accuracy).
% Junpeng Lao, 2015 August, University of Fribourg
%--------------------------------------------------------------------------
% Copyright (C) iMap Team 2015
%% Generate dataset - GMM
clear all;clc;
p = userpath;
addpath(genpath([ p(1:end-1) '/Apps/iMAP']));

nsamp = 10000;
% save current random generation state
defaultStream = RandStream.getGlobalStream();
savedState    = defaultStream.State;
% use the Multiplicative Lagged Fibonacci algotrithm for independent substreams
mystream      = RandStream.create('mlfg6331_64','NumStreams',nsamp,'StreamIndices',1);

RandStream.setGlobalStream(mystream);
reset(defaultStream); % allows restarting allays the same

xSize = 100;
ySize = 100;

muG          = [40,55; 60,55; 50,30; 50,48];
sigma(:,:,1) = [35 0; 0 30];
sigma(:,:,2) = [35 0; 0 30];
sigma(:,:,3) = [50 0; 0 60];
sigma(:,:,4) = [35 0; 0 60];
p            = [30,30,50,10];% the mixing parameter for GMM
obj          = gmdistribution(muG,sigma,p);

subplot(1,3,1)
ezsurf(@(x,y)pdf(obj,[x y]),[0 xSize],[0 ySize])
zlim([0,10])

subplot(1,3,2)
Nfix = 10; % total number of fixation
Y    = random(obj,Nfix);

plot(Y(:,1),Y(:,2),'o')
axis([0 xSize 0 ySize],'square')
%
% smooth map
smoothingpic   = 5;
[x, y]         = meshgrid(-floor(ySize/2)+.5:floor(ySize/2)-.5, -floor(xSize/2)+.5:floor(xSize/2)-.5);
gaussienne     = exp(- (x .^2 / smoothingpic ^2) - (y .^2 / smoothingpic ^2));
gaussienne     = (gaussienne - min(gaussienne(:))) / (max(gaussienne(:)) - min(gaussienne(:)));
% f_fil          = fft2(gaussienne);
% fixation matrix
coordX         = round(Y(:,2));
coordY         = round(Y(:,1));
intv           = normrnd(0.4,.085,length(Y),1);
indx1          = coordX>0 & coordY>0 & coordX<xSize & coordY<ySize;
rawmap         = full(sparse(coordX(indx1),coordY(indx1),intv(indx1),ySize,xSize));

% f_mat          = fft2(rawmap); % 2D fourrier transform on the points matrix
% filtered_mat   = f_mat .* f_fil;
% smoothpic      = real(fftshift(ifft2(filtered_mat)));
smoothpic      = conv2(rawmap, gaussienne,'same');

subplot(1,3,3)
imagesc(smoothpic);colorbar
set(gca,'YDir','normal');
axis('square','off')
%% Dataset generation
Ns          = 10;
Group       = {'CN','PA'};% control group and patient group
Ntrial      = 25;
MeanNfix    = 14;
stdNfix     = 3;
Meandur     = .4;
Stddur      = .085;
MC          = 0;
itt         = 0;
descriptemp = zeros(Ns*length(Group)*Ntrial, 10);
FixMap      = zeros(Ns*length(Group)*Ntrial, ySize, xSize);
RawMap      = FixMap;

for ig = 1:length(Group)
    %%
    figure;
    vidObj = VideoWriter(char(['Group' num2str(ig) '.avi']));
    open(vidObj)
    for is = 1:Ns
        % set seed
        MC       = MC+1;
        mystream = RandStream.create('mlfg6331_64','NumStreams',nsamp,'StreamIndices',MC);
        RandStream.setGlobalStream(mystream);
        if ig == 1
            % for Control group, only one type of generative model
            p      = [30,30,50,10];% the mixing parameter for GMM
            shif   = 5-randi(10,1,4);
            shifmn = 5-randi(10,size(muG));
            shifsd = randi(10,size(sigma));shifsd(2,1,:)=shifsd(1,2,:);
            pnew   = p+shif;
            obj    = gmdistribution(muG+shifmn,sigma+shifsd,pnew);
        else
            % for patient group, two generative model
            p1     = [30,30,50,10];% the mixing parameter for GMM
            shif   = 5-randi(10,1,4);
            shifmn = 5-randi(10,size(muG));
            shifsd = randi(10,size(sigma));shifsd(2,1,:)=shifsd(1,2,:);
            pnew1  = p1+shif;
            obj1   = gmdistribution(muG+shifmn,sigma+shifsd,pnew1);
            p2     = [6,6,50,10];% the mixing parameter for GMM
            shif   = 5-randi(10,1,4);
            shifmn = 5-randi(10,size(muG));
            shifsd = randi(10,size(sigma));shifsd(2,1,:)=shifsd(1,2,:);
            pnew2  = p2+shif;
            obj2   = gmdistribution(muG+shifmn,sigma+shifsd,pnew2);
        end
        if ig == 1
            ACC = (rand(1,Ntrial)>.6)+1;% 1 correct, 2 incorrect
        else
            ACC = (rand(1,Ntrial)>.4)+1;% 1 correct, 2 incorrect
        end
        % ACC=(rand(1,Ntrial)>.5)+1;% 1 correct, 2 incorrect
        ACCtmp  = rand(size(ACC));
        [a,b]   = sort(ACC);
        ACC2    = zeros(size(ACC));
        ACC2(b) = 1-sort(ACCtmp);
        for it = 1:Ntrial
            itt  = itt+1;
            Nfix = ceil(normrnd(MeanNfix,stdNfix)); % total number of fixation
            if ig == 2
                if ACC(it) == 1;
                    obj = obj1;
                else
                    obj  = obj2;
                    Nfix = ceil(normrnd(MeanNfix*.78,stdNfix)); % total number of fixation
                end
            end
            
            Ytmp   = random(obj,Nfix);
            Ytmp2  = [randi(xSize,2,1) randi(ySize,2,1)];
            Y      = [Ytmp;Ytmp2];
            hold on
            plot(Y(:,1),Y(:,2),'.','color',[0 0 0])
            drawnow
            axis([0 xSize 0 ySize],'square','off')
            currFrame = getframe;
            writeVideo(vidObj,currFrame);
            
            coordX     = xSize-round(Y(:,2));
            coordY     = round(Y(:,1));
            pathlength = diag(squareform(pdist([coordY,coordX])),1);
            intv       = normrnd(Meandur,Stddur,length(Y),1)*1000;
            indx1      = coordX>0 & coordY>0 & coordX<xSize & coordY<ySize;
            rawmap     = full(sparse(coordX(indx1),coordY(indx1),intv(indx1),ySize,xSize));
            
%             f_mat              = fft2(rawmap); % 2D fourrier transform on the points matrix
%             filtered_mat       = f_mat .* f_fil;
%             smoothpic          = real(fftshift(ifft2(filtered_mat)));
            smoothpic          = conv2(rawmap, gaussienne,'same');
            mm                 = mean(smoothpic(:));
            stdm               = std(smoothpic(:));
            FixMap(itt,:,:)    = (smoothpic-mm)./stdm;
            RawMap(itt,:,:)    = rawmap;
            descriptemp(itt,:) = [Nfix, ...
                sum(intv), sum(intv)/Nfix, sum(pathlength),mean(pathlength), ...
                ig,        ACC(it),        is+Ns*(ig-1),   it, ACC2(it)];
        end
    end
    close(vidObj);
end
%%
table_header2     = [{'FixNum'},{'sumFixDur'},{'meanFixDur'},{'totalPathLength'},...
    {'meanPathLength'},{'Grp'},{'ACC'},{'Sbj'},{'Trial'},{'ACC2'}];
DescriptvM1       = [table_header2; num2cell(descriptemp)];
DescriptvM1       = cell2dataset(DescriptvM1);

DescriptvM1.Trial = nominal(DescriptvM1.Trial);
DescriptvM1.Sbj   = nominal(DescriptvM1.Sbj);
DescriptvM1.Grp   = nominal(DescriptvM1.Grp, Group);
DescriptvM1.ACC   = nominal(DescriptvM1.ACC, {'hit','miss'});

PredictorM        = DescriptvM1(:,6:end);
DescriptvM        = DescriptvM1(:,1:end-2);

Mask              = squeeze(mean(FixMap,1))>.1;
%% save matrix
save(strcat('./FixMap_single_trial_scaled'), 'FixMap',     '-v7.3');
save(strcat('./PredictorM_single_trial'),    'PredictorM', '-v7.3');
save(strcat('./DescriptvM_single_trial'),    'DescriptvM', '-v7.3');
save(strcat('./RawMap_single_trial_scaled'), 'RawMap',     '-v7.3');
save(strcat('./Mask_single_trial_scaled'),   'Mask',       '-v7.3');
%% mean map
descriptive_part(DescriptvM, FixMap, 0)
CondiVec     = PredictorM.Grp .* PredictorM.ACC;
SbjVec       = PredictorM.Sbj;
[RDM, stRDM] = rdmfixmap(FixMap, Mask, CondiVec, SbjVec);
%% LMM
tic
opt.singlepredi    = 1;
% opt.parallelname='grid2';
[LMMmap,lmexample] = imapLMM(FixMap,PredictorM,Mask,opt, ...
    'PixelIntensity ~ Grp * ACC + (1|Sbj)', ...
    'DummyVarCoding','effect');
save('LMMmap_ACC.mat','LMMmap','-v7.3');
toc
%% plot model fitting
opt1      = struct; 
opt1.type = 'model';
% perform contrast
[StatMap] = imapLMMcontrast(LMMmap,opt1);
% output figure;
imapLMMdisplay(StatMap,0)
%% plot fixed effect (anova result using the cell mean DS and its related contrast)
% close all
opt        = struct;% clear structure
opt.type   = 'predictor beta';
opt.c      = limo_OrthogContrasts([2,2]);
opt.name   = {'Grp','ACC','Interaction'};
opt.alpha  = 0.05;
% perform contrast
[StatMap]  = imapLMMcontrast(LMMmap,opt);
imapLMMdisplay(StatMap,0);

mccopt           = struct;
mccopt.methods   = 'bootstrap';
mccopt.bootopt   = 1;
mccopt.bootgroup = {'Grp'};
mccopt.nboot     = 1000;
[StatMap_c]      = imapLMMmcc(StatMap,LMMmap,mccopt,FixMap);
imapLMMdisplay(StatMap_c,0);

%% post-hoc
[Posthoc]        = imapLMMposthoc(StatMap_c,RawMap,LMMmap,'mean')
%% LMM 2
tic
opt                 = struct;% clear structure
opt.singlepredi     = 1;
% opt.parallelname='grid2';
[LMMmap2,lmexample] = imapLMM(FixMap,PredictorM,Mask,opt, ...
    'PixelIntensity ~ Grp * ACC2 + (1|Sbj)', ...
    'DummyVarCoding','effect');
save('LMMmap_ACC2.mat','LMMmap2','-v7.3');
toc
%% plot model fitting
opt1      = struct;% clear structure
opt1.type = 'model';
% perform contrast
[StatMap] = imapLMMcontrast(LMMmap2,opt1);
% output figure;
imapLMMdisplay(StatMap,0)
%% plot fixed effect (anova result using the cell mean DS and its related contrast)
close all
opt       = struct;% clear structure
opt.type  = 'fixed';
% perform contrast
[StatMap] = imapLMMcontrast(LMMmap2,opt);
imapLMMdisplay(StatMap,0);

mccopt           = struct;
mccopt.methods   = 'bootstrap';
mccopt.bootopt   = 1;
mccopt.bootgroup = {'Grp'};
mccopt.nboot     = 1000;
[StatMap_c]      = imapLMMmcc(StatMap,LMMmap2,mccopt,FixMap);
imapLMMdisplay(StatMap_c,0);

%%
opt       = struct;% clear structure
opt.type  = 'model beta';
% perform contrast
[StatMap] = imapLMMcontrast(LMMmap2,opt);
imapLMMdisplay(StatMap,0);

mccopt           = struct; 
mccopt.methods   = 'bootstrap';
mccopt.bootopt   = 1;
mccopt.bootgroup = {'Grp'};
mccopt.nboot     = 1000;
[StatMap_c]      = imapLMMmcc(StatMap,LMMmap2,mccopt,FixMap);
imapLMMdisplay(StatMap_c,0);
