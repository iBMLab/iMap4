% the following script is a demonstration of the Single-Trial estimation
% method for iMap4. 
% A computer simulated data set is introduced here: a 4 x 4 grid was
% presented to 20 subjects for 100 trials, each trial participant gives a
% subjective rating (from -3 to 3 in the current case). We simulate a
% linear relationship between rating and fixation number on each grid with
% different power. iMap4 will then estimate the relationship for each
% grid.
% Junpeng Lao, 2015 March, University of Fribourg
%--------------------------------------------------------------------------
% Copyright (C) iMap Team 2015
%% create linear relationship sample - one subject
p = userpath;
addpath(genpath([ p(1:end-1) '/Apps/iMAP']));
rng('default'); % For reproducibility
rng(1)

% introduce linear relationship 
rho     = [0.9, 0.6,  0.3,  0  ];
slope   = [1.0, 0.4, -0.2, -0.8];
n       = 100; % N trials
iplot   = 0;

mu1     = [0 0];
Sigma1  = [1 0; 0 1];
R1      = chol(Sigma1);
Z1      = repmat(mu1,n,1) + randn(n,2)*R1;
X       = Z1(:,1);

figure;

for irr = 1:length(rho)
    for iss = 1:length(slope)
        % Gaussian copula
        % mu       = [0 0];
        Sigma    = [1 rho(irr); rho(irr) 1];
        % Z        = mvnrnd(mu,Sigma,n);
        % U        = normcdf(Z,0,1);
        R        = chol(Sigma);
        Z1(:,2)  = randn(n,1);
        U        = Z1*R;
        
        iplot    = iplot+1;
        subplot(4,4,iplot);hold on
        % introduce slope
        U(:,2)   = U(:,2).*slope(iss);
        % normalized
        constant = min(U(:,2));
        cst      = constant * sign(constant);
        U(:,2)   = U(:,2) + cst;
        
        plot(X,U(:,2),'o');
        plot(X,X.*slope(iss)+cst,'r','LineWidth',2.5)
        title(['Y_{' num2str(iplot) '}=' num2str(slope(iss)) '*X+c_{' num2str(iplot) '}; (r_{' num2str(iplot) '}=' num2str(rho(irr)) ')']);
        axis([-2.2 2.2 0 5.5])
    end
end
xlabel('Rating');
ylabel('Response');
%% create data sample - all subject
figure;

Ns        = 20;% Ns subjects
mu1       = [0 0];
Sigma1    = [1 0; 0 1];
R1        = chol(Sigma1);
Xall      = zeros(Ns, 100);
Yall      = zeros(Ns, length(rho)*length(slope), 100);
for is = 1:Ns
    rng(is)
    Z1         = repmat(mu1,n,1) + randn(n,2)*R1;
    Xall(is,:) = Z1(:, 1);
    itype      = 0;
    for irr = 1:length(rho)
        for iss = 1:length(slope)
            itype    = itype+1;
            % Gaussian copula
            Sigma    = [1 rho(irr); rho(irr) 1];
            R        = chol(Sigma);
            Z1(:,2)  =  randn(n,1);
            U        = Z1*R;
            
            % introduce slope
            U(:,2)   = U(:,2).*slope(iss);
            % normalized
            constant = min(U(:,2));
            cst      = constant*sign(constant);
            U(:,2)   = U(:,2)+cst;
            Yall(is,itype,:) = U(:,2);
            if itype==1
                plot(Xall(is,:),squeeze(Yall(is,1,:)),'o')
                hold on
            end
        end
    end
end
%% generate fixation matrix according to the design matrix using a Gaussian mixture model
figure;
xSize    = 150;
ySize    = 150;

ip = 0;
for iy = 1:length(rho)
    for ix = 1:length(slope)
        ip = ip+1;
        muG(ip,:) = [ix*30,iy*30];
    end
end
sigma        = [10 0; 0 10];
p            = squeeze(Yall(1,:,7));% the mixing parameter for GMM
Nfix         = 60;
muG2         = muG;
muG2(p==0,:) = [];
p   (p==0)   = [];

obj          = gmdistribution(muG2,sigma,p);
subplot(1,3,1)
ezsurf(@(x,y)pdf(obj,[x y]),[0 xSize],[0 ySize])
zlim([0,10])

subplot(1,3,2)
rng(1); % For reproducibility
Nfix         = 60; % total number of fixation
Y            = random(obj,Nfix);

plot(Y(:,1),Y(:,2),'o')
set(gca,'YDir','reverse');
axis([0 xSize 0 ySize],'square')
% smooth map
smoothingpic = 5;
[x, y]       = meshgrid(-floor(ySize/2)+.5:floor(ySize/2)-.5, -floor(xSize/2)+.5:floor(xSize/2)-.5);
gaussienne   = exp(- (x .^2 / smoothingpic ^2) - (y .^2 / smoothingpic ^2));
gaussienne   = (gaussienne - min(gaussienne(:))) / (max(gaussienne(:)) - min(gaussienne(:)));
f_fil        = fft2(gaussienne);
% fixation matrix
coordX       = round(Y(:,2));
coordY       = round(Y(:,1));
indx1        = coordX>0 & coordY>0 & coordX<xSize & coordY<ySize;
rawmap       = full(sparse(coordX(indx1),coordY(indx1),ones(size(coordY(indx1))),ySize,xSize));
f_mat        = fft2(rawmap); % 2D fourrier transform on the points matrix
filtered_mat = f_mat .* f_fil;
smoothpic    = real(fftshift(ifft2(filtered_mat)));
subplot(1,3,3)
imagesc(smoothpic);colorbar
axis('square','off')

%% create FixMap and PredictorM
Subject      = zeros(Ns*n,1);
Rating       = zeros(Ns*n,1);
fixNb        = zeros(Ns*n,1);
FixMap       = zeros(Ns*n,ySize,xSize);
RawMap       = FixMap;
Nfixbase     = 2; % number of fixations for each trials is the multiple of Nfixbase
item         = 0;
DataM        = [];
for is = 1:Ns
    for itrial = 1:n
        item = item+1;
        rng(item)
        p            = squeeze(Yall(is,:,itrial));% the mixing parameter for GMM
        Nfix         = round(sum(p))*Nfixbase;
        muG2         = muG;
        muG2(p==0,:) = [];
        p   (p==0)   = [];
        obj          = gmdistribution(muG2,sigma,p);
        Ytmp         = random(obj,Nfix);
        Ytmp2        = [randi(xSize,10,1) randi(ySize,10,1)];% 10 random fixations 
        % fixation matrix
        Y            = [Ytmp;Ytmp2];
        coordX       = round(Y(:,2));
        coordY       = round(Y(:,1));
        indx1        = coordX>0 & coordY>0 & coordX<xSize & coordY<ySize;
        rawmap       = full(sparse(coordX(indx1),coordY(indx1),ones(size(coordY(indx1))),ySize,xSize));
        f_mat        = fft2(rawmap); % 2D fourrier transform on the points matrix
        filtered_mat = f_mat .* f_fil;
        smoothpic    = real(fftshift(ifft2(filtered_mat)));
        
        % save fixation matrix
        RawMap (item, :, :) = rawmap;
        FixMap (item, :, :) = smoothpic;
        Subject(item)       = is;
        Rating (item)       = Xall(is,itrial);
        fixNb  (item)       = Nfix;
        
        Yk                  = rawmap(rawmap~=0);        
        [Yk(:,2),Yk(:,3)]   = find(rawmap~=0);
        Yk(:,4)             = Xall(is,itrial);
        Yk(:,5)             = is;
        DataM               = [DataM; Yk];
    end
end
% figure;imagesc(squeeze(mean(FixMap,1)))
Mask        = squeeze(mean(FixMap,1)>.25);
% figure;imagesc(Mask)
PredictorM  = dataset(Rating,fixNb,Subject);

ds          = mat2dataset(DataM);
ds          = set(ds,'VarNames',{'fixN','rl','cl','rating','is'});
export(ds,'File','STrate.csv','Delimiter',',');
%% Linear Mixed Modeling with imapLMM
tic
opt.singlepredi    = 0;
[LMMmap,lmexample] = imapLMM(FixMap,PredictorM,Mask,opt, ...
    'PixelIntensity ~ Rating + (1|Subject)', ...
    'DummyVarCoding','effect','FitMethod','REML');
toc
%% show result
%% plot model fitting
close all
opt1.type = 'model';
% perform contrast
[StatMap] = imapLMMcontrast(LMMmap,opt1);
% output figure;
imapLMMdisplay(StatMap,0);
%% plot fixed effec(anova and beta)
% close all
opt       = struct;% clear structure
opt.type  = 'model beta';
opt.alpha = 0.05;
% perform contrast
[StatMap] = imapLMMcontrast(LMMmap,opt); 
% imapLMMdisplay(StatMap,1)

mccopt         = struct;
mccopt.methods = 'bootstrap';
mccopt.nboot   = 1000;
mccopt.bootopt = 1;
mccopt.tfce    = 0;
% perform multiple comparison correction
[StatMap_c]    = imapLMMmcc(StatMap,LMMmap,mccopt,FixMap);
% output figure;
imapLMMdisplay(StatMap_c,0)
%% Region dependent threshold
ResampStat     = StatMap_c.resampMat;
nboot          = mccopt.nboot;
alpha          = 0.05;
Pmask          = zeros(size(squeeze(StatMap_c.Pmask(2,:,:))));
Pmap           = StatMap.Pmap;
mapvalue       = StatMap.map;

Maskall        = zeros(size(Pmask));
Maskall(16:135,16:135) = kron(reshape(1:16,4,4),ones(30,30));
bootopt = 1;

for imask = 1:16
    region = Maskall==imask;
    for ic = 2
        Fboot  = squeeze(ResampStat.resFvalue(:,ic,region));
        Pboot  = squeeze(ResampStat.resPvalue(:,ic,region));
        Betabt = squeeze(ResampStat.resBeta  (:,ic,region));
        Fboot  = reshape(Fboot,  1000, 30, 30);
        Pboot  = reshape(Pboot,  1000, 30, 30);
        Betabt = reshape(Betabt, 1000, 30, 30);
        % output estimated bootstrap distribution (cluster test)
        maxclust = zeros(nboot,3);
        for ib = 1:nboot
            Ftmp            = squeeze(Fboot(ib,:,:));
            [bmasktmp,bnum] = bwlabel(squeeze(Pboot(ib,:,:))<alpha);
            if bnum > 0
                maxtmp = zeros(bnum, 3);
                for icluster = 1:bnum
                    maxtmp(icluster, 1) = nansum(nansum(Ftmp(bmasktmp==icluster)));
                    maxtmp(icluster, 2) = sum(sum(bmasktmp==icluster));
                    maxtmp(icluster, 3) = maxtmp(icluster,1)./maxtmp(icluster,2);
                end
                maxclust(ib,1) = max(maxtmp(:,1));
                maxclust(ib,2) = max(maxtmp(:,2));
                maxclust(ib,3) = max(maxtmp(:,3));
            end
        end
        distmp  = sort(maxclust,1);
        % output cluster mass and size threshold at alpha
        clthres = distmp(round(nboot*(1-alpha)),:);
        % output new Pmask: 1 cluster mass, 2 cluster size, 3 both,
        % 4 cluster dense
        mapvalue1 = reshape(squeeze(mapvalue(ic,region)), 30, 30);
        Pmap1     = reshape(squeeze(Pmap(ic,region)),     30, 30);
        switch bootopt
            case 1
                Pmask(region) = clustertest2D(mapvalue1, ...
                    Pmap1,alpha,clthres(1),[],[]);
            case 2
                Pmask(region) = clustertest2D(mapvalue1, ...
                    Pmap1,alpha,[],clthres(2),[]);
            case 3
                Pmask(region) = clustertest2D(mapvalue1, ...
                    Pmap1,alpha,clthres(1),clthres(2),[]);
            case 4
                Pmask(region) = clustertest2D(mapvalue1, ...
                    Pmap1,alpha,[],[],clthres(3));
        end
    end
end
figure;imagesc(Pmask);axis square off
%% plot fixed effec(anova and beta) - using permutation
% close all
opt       = struct;% clear structure
opt.type  = 'model beta';
opt.alpha = .05;
% perform contrast
[StatMap] = imapLMMcontrast(LMMmap,opt); 
% imapLMMdisplay(StatMap,1)

mccopt         = struct;
mccopt.methods = 'permutation';
mccopt.nboot   = 1000;
mccopt.bootopt = 1;
mccopt.tfce    = 0;
% perform multiple comparison correction
[StatMap_c]    = imapLMMmcc(StatMap,LMMmap,mccopt,FixMap);
% output figure;
imapLMMdisplay(StatMap_c,0)
%% Region dependent threshold
ResampStat     = StatMap_c.resampMat;
nboot          = mccopt.nboot;
alpha          = 0.05;
Pmask          = zeros(size(squeeze(StatMap_c.Pmask(2,:,:))));
Pmap           = StatMap.Pmap;
mapvalue       = StatMap.map;

Maskall        = zeros(size(Pmask));
Maskall(16:135,16:135) = kron(reshape(1:16,4,4),ones(30,30));
bootopt  = 1;
Pmapnew2 = Maskall;
for imask = 1:16
    region = Maskall==imask;
    for ic = 2
        Fboot    = squeeze(ResampStat.resFvalue(:,ic,region));
        Fboot    = reshape(Fboot,1000,30,30);
        % output new Pmask
        origFmap         = squeeze(StatMap.map(ic,region));
        origFmap         = reshape(origFmap,30,30);
        orFmat           = permute(repmat(origFmap,[1,1,nboot]),[3,1,2]);
        Fboot2           = repmat(max(Fboot(:,:),[],2),[1,size(Fboot,2),size(Fboot,3)]);
        Pmapnew2(region) = (sum(Fboot2>=orFmat,1)+1)./(nboot+1);
    end
end
Pmapnew2( Mask==0 ) = 1;
figure;
subplot(1,3,1);
imagesc(squeeze(StatMap_c.Pmap(ic,:,:)));
axis square off;
subplot(1,3,2);
imagesc(Pmapnew2);
axis square off;
subplot(1,3,3);
imagesc(Pmapnew2<alpha);
axis square off;
