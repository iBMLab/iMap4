%--------------------------------------------------------------------------
% code snippet to demonstrate the theory of linear mixed model and
% permutation for LMM.
% Copyright (C) iMap Team 2015
%% Permutation on LMM
% parameters
Ns       = 10; % number of subject
Condi1   = 2;  % levels of condition 1
Condi2   = 3;  % levels of condition 2
Nt       = 10; % number of trials
beta0    = [1 2 3,2 2 1]';% [g1_c1 g1_c2 g1_c3, g2_c1 g2_c2 g2_c3]
sbjMean  = 0;  
sbjVar   = 1.5;
noisevar = 5;
% Generate dataset
DXss         = kron(eye(6),ones(Nt,1));
DXall        = repmat(DXss,[Ns,1]);
% sbjidx=reshape(repmat(1:Ns,[Nt*Condi1*Condi2,1]),[Ns*Nt*Condi1*Condi2,1]);
sbjidx       = reshape(repmat(1:(Ns*2),[Nt*Condi2,1]),[Ns*Nt*Condi1*Condi2,1]);
sbjintercept = random('Normal',sbjMean,sbjVar,Ns,1);
sbjintM      = zeros(size(sbjidx));
for is=1:Ns
    sbjintM(sbjidx==is)=sbjintercept(is);
end
% Y=random('Normal',0,noisevar,[Ns*Condi1*Condi2*Nt,1])+DXall*beta0+sbjintM;
Y = random('gp',1,1,0,[Ns*Condi1*Condi2*Nt,1])./6+DXall*beta0+sbjintM;
%
% HLM
beta1 = zeros(Ns,Condi2);
DXss2 = DXss(1:Nt*Condi2,1:Condi2);
for is=1:Ns
    Ytmp        = Y(sbjidx==is);
    DXtmp       = [DXss2,ones(length(DXss2),1)];
    % DXtmp=DXss2;
    invD        = pinv(DXtmp);
    betatmp     = invD*Ytmp;
    beta1(is,:) = betatmp(1:Condi2);
end
DX2tmp = repmat(1:Condi2,Ns,1);DX2tmp(2:2:Ns,:)=DX2tmp(2:2:Ns,:)+Condi2;
Y2     = beta1(:);DX2=DX2tmp(:);
DXg    = dummyvar(DX2);
dfehlm = length(Y2)-rank(DXg);
beta2  = DXg\Y2;
mse    = sum((DXg*beta2-Y2).^2)/dfehlm;
covhlm = mse*((DXg'*DXg)^-1);

%
% LMM 2
lmm2    = fitlmematrix(DXall,Y,ones(length(Y),1),sbjidx,'Dummyvarcoding','effect','FitMethod', 'REML');%,'CovariancePattern','Diagonal');
betalmm = lmm2.Coefficients.Estimate;
covlmm  = lmm2.CoefficientCovariance;

% LMM 1
tbl        = dataset;
tbl.Y      = Y;
tbl.Condi1 = nominal(2-double(sum(DXall(:,1:3),2)>0));
tbl.Condi2 = nominal(double(sum(DXall(:,[1 4]),2)>0)+double(sum(DXall(:,[2 5]),2)>0)*2+double(sum(DXall(:,[3 6]),2)>0)*3);
tbl.sbjidx = sbjidx;
lmm1       = fitlme(tbl,'Y~Condi1*Condi2+(1|sbjidx)','Dummyvarcoding','effect','FitMethod', 'REML');
anova(lmm1)

barall = zeros(2,3);
for ic1=1:2
    for ic2=1:3
        indx = double(tbl.Condi1)==ic1&double(tbl.Condi2)==ic2;
        barall(ic1,ic2) = mean(Y(indx));
    end
end
% plot result
figure;
subplot(1,3,1)
bar(mean(barall,2))

subplot(1,3,2)
bar(mean(barall,1))

subplot(1,3,3)
bar(barall)

% ANOVA
C = limo_OrthogContrasts([Condi1,Condi2]);
clear pl Fl df1l df2l ph Fh df1h df2h
for ic=1:3
    c = C{ic};
    [pl(ic),Fl(ic),df1l(ic),df2l(ic)] = coefTest(lmm2,c);
    df1h(ic) = rank(c);
    df2h(ic) = dfehlm;%(Ns-1)*rank(c);
    Fh(ic)   = ((c*beta2)'*((c*covhlm*c')^-1)*(c*beta2))./rank(c);
    ph(ic)   = 1-fcdf(Fh(ic),df1h(ic),df2h(ic));
end
FtableLME = dataset(Fl',df1l',df2l',pl','Varnames',{'Fstat','DF1','DF2','pValue'},'ObsNames',{'Condi1','Condi2','Interation'})
FtableHLM = dataset(Fh',df1h',df2h',ph','Varnames',{'Fstat','DF1','DF2','pValue'},'ObsNames',{'Condi1','Condi2','Interation'})

%% Matrix partitioning for interaction
DX1=lmm1.designMatrix; % Sigma-restricted Designmatrix
DX2=lmm2.designMatrix; % cell mean Designmatrix
mdlopt=2;
figure;
subplot(2,3,1);imagesc(DX1)
subplot(2,3,2);imagesc(DX2)

Nperm=10000;
if mdlopt==1
    DX   = DX1;
    Y    = lmm1.response;
    Zx   = designMatrix(lmm1,'Random');
    BLUP = randomEffects(lmm1);
    bY   = Y - Zx*BLUP;
    fitml= strcmp(lmm1.FitMethod,'ML');
    % Partition matrix for testing interaction and output original stats
    c2   = [0 0 0 0 1 0;0 0 0 0 0 1];
    [p,F,df1,df2]=coefTest(lmm1,c2);
    disp(['F test for the interaction (SR-DX): F(',num2str(df1),',',num2str(df2) ') = ', num2str(F), '; p = ', num2str(p)])
else
    DX   = DX2;
    Y    = lmm2.response;
    Zx   = designMatrix(lmm2,'Random');
    BLUP = randomEffects(lmm2);
    bY   = Y - Zx*BLUP;
    fitml= strcmp(lmm1.FitMethod,'ML');
    % Partition matrix for testing interaction and output original stats
    c2   = C{ic};
    [p,F,df1,df2]=coefTest(lmm2,c2);
    disp(['F test for the interaction (MC-DX): F(',num2str(df1),',',num2str(df2) ') = ', num2str(F), '; p = ', num2str(p)])
end


% the following part is from Appendix A. of Winkler, et al (2014).
% Permutation inference for the general linear model. Neuroimage, 92,
% 381-397.

%%%%%%%%%% This part of the code is adapted from imapLMMresample %%%%%%%%%%
% Partition design matrix
cu     = null(c2,'r');
c2     = c2';
D      = pinv(DX'*DX);
cv     = cu-c2*pinv(c2'*D*c2)*c2'*D*cu;
parX   = DX*D*c2*pinv(c2'*D*c2);
parZ   = DX*D*cv*pinv(cv'*D*cv);
M      = [parX,parZ]; % new desgin matrix
cnew   = zeros(size(c2'));
for ir=1:rank(parX);cnew(ir,ir)=1;end
% orignial statistic for the fixed effect (revisit)
% bY2=(1-parZ*pinv(parZ))*bY;
obeta  = M\bY;
cz     = ones(1,rank(M))-sum(cnew,1);
bY2    = bY-M(:,cz==1)*obeta(cz==1,:);
% remove nuisance effects
bbeta  = M\bY2;
if fitml
    bMSE   = mean((M*bbeta-bY2).^2); % ML
else
    bMSE   = sum((M*bbeta-bY2).^2)./df2; % REML
end
covlmm = bMSE*((M'*M)^-1);
Fparti1= ((cnew*bbeta)'*((cnew*covlmm*cnew')^-1)*(cnew*bbeta))./df1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% This part of the code is adapted from MRM_estimateContrastsPerm %%%%%
% Model partition after Ridgway (2009) - see Winkler et al. (2014)
X     = DX * pinv(c2');
Z     = DX - (DX * c2 * pinv(c2));
[Z,~] = svd(Z);
Z     = Z(:,1:rank(DX) - rank(c2));
Rz    = eye(size(Z,1)) - Z * pinv(Z);    % (I - ZZ*)
X     = Rz * X;                          % Orthogonalise X wrt Z
M2    = [X Z];
% New contrast for testing K = 0 in Y* = XK + ZG + E*, specified in
% terms of the combined model Y* = [X Z]D = E*
Con   = padarray(eye(size(c2',1)), [0 (size(X,2) + size(Z,2)) - size(c2',1)], 0, 'post');

% Notice, the above is the same compare to the function written by Winkler
% [X,Z,eCm,eCx] = palm_partition(DX,c2,'ridgway');
% M2    = [X Z];
% Con   = eCm';

pXX   = pinv(M2' * M2);            % (M'M)*
CXXC  = pinv(Con * pXX * Con');    % (C (M'M)* C')*
pM    = pinv(M2);                  %  M*
Ry    = eye(size(M2,1)) - M2 * pM; % (I - MM*)
ResidZ = Rz * bY;
PEs    = pM * ResidZ;
Resids = Ry * ResidZ;
bhat   = Con*PEs;
if fitml
    Fparti2= ((bhat'*(X'*X)*bhat)/rank(Con))/((Resids'*Resids)/(size(M2,1)));
else
    Fparti2= ((bhat'*(X'*X)*bhat)/rank(Con))/((Resids'*Resids)/(size(M2,1)-rank(X)-rank(Z)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% permutation
Fperm1=zeros(Nperm,1);
Fperm2=zeros(Nperm,1);
tic
srand = rng;
for iperm=1:Nperm
    permlist = randperm(size(M,1),size(M,1));
    
    bbeta  = M\bY2(permlist);
    if fitml
        bMSE   = mean((M*bbeta-bY2(permlist)).^2);
    else
        bMSE   = sum((M*bbeta-bY2(permlist)).^2)./df2;
    end
    covlmm = bMSE*((M'*M)^-1);
    Fperm1(iperm) = ((cnew*bbeta)'*((cnew*covlmm*cnew')^-1)*(cnew*bbeta))./df1;
end
time1 = toc;
p1=mean(Fperm1>=Fparti1);
disp(['Computational time for ' num2str(Nperm) ' permutation is ' num2str(time1)])
disp(['Permutation F test 1 for the interaction: F(',num2str(df1),',',num2str(df2) ') = ', num2str(Fparti1), '; p = ', num2str(p1)])

tic
rng(srand);
for iperm=1:Nperm
    permlist = randperm(size(M,1),size(M,1));

    PEs    = pM * ResidZ(permlist);
    Resids = Ry * ResidZ(permlist);
    bhat   = Con*PEs;
    if fitml
        Fperm2(iperm) = ((bhat'*(X'*X)*bhat)/rank(Con))/((Resids'*Resids)/(size(M2,1)));
    else
        Fperm2(iperm) = ((bhat'*(X'*X)*bhat)/rank(Con))/((Resids'*Resids)/(size(M2,1)-rank(X)-rank(Z)));
    end
end
time2 = toc;
p2=mean(Fperm2>=Fparti2);
disp(['Computational time for ' num2str(Nperm) ' permutation is ' num2str(time2)])
disp(['Permutation F test 2 for the interaction: F(',num2str(df1),',',num2str(df2) ') = ', num2str(Fparti2), '; p = ', num2str(p2)])

subplot(2,3,4);imagesc(DX);colorbar
subplot(2,3,5);imagesc(M);colorbar
subplot(2,3,6);imagesc(M2);colorbar

