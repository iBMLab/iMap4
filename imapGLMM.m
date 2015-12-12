function [GLMMmap,glmexample]=imapGLMM(FixMap,PredictorM,Mask,opt,formula,varargin)
% Usage: [GLMMmap,glmexample]=imapGLMM(FixMap,PredictorM,Mask,opt,formula,varargin)
% Input:    FixMap - total number of trials * xSize * ySize
%       PredictorM - dataset format of condition Matrix, total number of
%                    trials * number of predictor. Categorical column must
%                    set to nominal
%             Mask - 2D mask, for reducing the number of computation
%              opt - structure. option to define parallel grid (opt.parallelname)
%                    and option to compute each single categorical condition
%                    beta (opt.singlepredi)
% formula/varargin - same as fitglme, type '>>help fitglme' for more
%                    information
%
% Only apply when the input Predictor Matrix has non-orthogonal categorical
% predictor. Categorical beta of each condition will be compute separately
% when opt.singlepredi is set to 1.
% PredictorMatrix should include one subject colume to indicate group for
% the estimation of random effect.
% The function will run in Parallel if Parallel Computing Toolbox is
% installed in Matlab. Parallel cluster could be specify in
% opt.parallelname (default 'local')
% Same varargin apply as fitlme in matlab. For more information
% help fitlme
% See also fitlme, LinearMixedModel, LinearModel, GeneralizedLinearModel, NonLinearModel.
%
% 2014-11-08 Code updated to sort out compatible issue with Matlab 2014b
% 2015-02-12 Junpeng Lao, University of Fribourg.
%--------------------------------------------------------------------------
% Copyright (C) iMap Team 2015
%% check Matlab version and toolbox before start
v=ver;
if size(PredictorM,1)~=size(FixMap,1)
    error('Number of items/trials mismatch between PredictorM and FixMap, please double check your input matrix.')
end
if exist('GeneralizedLinearMixedModel')==0
    error('Some crucial toolboxes are missing in your current version of Matlab for imapLMM to run.')
end
has_fsolve=any(strcmp({v.Name},'Parallel Computing Toolbox'));
if has_fsolve==1
    if isfield(opt,'parallelname')==0
        gridname='local';
    else
        gridname=opt.parallelname;
    end
    disp('imapLMM will be computed in parallel mode, you machine might become quite slow during the computation.')
end
if isfield(opt,'singlepredi')==1;
    singlep=opt.singlepredi;
    if singlep~=0
        disp('Single predictor beta for each categorical condition will be computed.')
    end
end
% create save structure
GLMMmap=struct;
GLMMmap.runopt=opt;
% check mask
ySize2=size(FixMap,2);
xSize2=size(FixMap,3);
if isempty(Mask)==1
    Mask=ones(ySize2,xSize2);
end
compuIdx=find(Mask(:)==1);
numCompute=length(compuIdx);
% run one model as example to create matrix for saving
tic;
% indextmp=compuIdx(randi(numCompute));
valall=squeeze(mean(FixMap,1));
indextmp= valall(compuIdx)==max(valall(compuIdx));
tbl = PredictorM;% construct lme table
tbl.PixelIntensity = FixMap(:,indextmp);
if isa(tbl,'dataset')
    VarNames=tbl.Properties.VarNames;
elseif isa(tbl,'table')
    VarNames=tbl.Properties.VariableNames;
    tbl= table2dataset(tbl);
else
    error('Please reform your input predictor matrix as dataset!')
end

glmexample = GeneralizedLinearMixedModel.fit(tbl,formula,varargin{:});
[~,~,randstat]=randomEffects(glmexample);
time=toc;
disp(glmexample)
estimatime=time*numCompute/60;
disp(['The total computation time is around ' num2str(estimatime) ...
    ' minutes without parallel computation.'])
% save information
GLMMmap.Distribution=glmexample.Distribution;
GLMMmap.Link=glmexample.Link;

GLMMmap.VariableInfo=glmexample.VariableInfo;
GLMMmap.Variables=tbl;
GLMMmap.FitMethod=glmexample.FitMethod;
GLMMmap.Formula=glmexample.Formula;
GLMMmap.modelX=designMatrix(glmexample);
GLMMmap.FitOptions=varargin;
GLMMmap.modelDFE=glmexample.DFE;
GLMMmap.CoefficientNames=glmexample.CoefficientNames;
% save Anova information
statexample=anova(glmexample);% no need to save pValue=1-fcdf(FStat,DF1,DF2);
GLMMmap.Anova.Term=statexample.Term;
GLMMmap.Anova.DF1=statexample.DF1;
GLMMmap.Anova.DF2=statexample.DF2;

AnovaFStat=NaN(size(statexample.Term,1),ySize2,xSize2);
Covb=NaN([size(glmexample.CoefficientCovariance),ySize2,xSize2]);
MSEmatrix=NaN(ySize2,xSize2);
Dispersion=MSEmatrix;
SSEmatrix=MSEmatrix;
SSTmatrix=MSEmatrix;
SSRmatrix=MSEmatrix;
RsquaredMatrix=NaN(2,ySize2,xSize2);
ModelCriterionMatrix=NaN(4,ySize2,xSize2);% AIC,BIC,LogLikelihood,Deviance
CoefficientsMatrix=NaN(length(glmexample.CoefficientNames),2,ySize2,xSize2);% save beta, SE
% for parallel computing
pAnovaFStat=NaN(size(statexample.Term,1),numCompute);
pCovb=NaN([size(glmexample.CoefficientCovariance),numCompute]);
pMSEmatrix=NaN(1,numCompute);
pDispersion=NaN(1,numCompute);
pSSEmatrix=pMSEmatrix;
pSSTmatrix=pMSEmatrix;
pSSRmatrix=pMSEmatrix;
pRsquaredMatrix=NaN(2,numCompute);
pModelCriterionMatrix=NaN(4,numCompute);% AIC,BIC,LogLikelihood,Deviance
pCoefficientsMatrix=NaN(length(glmexample.CoefficientNames),2,numCompute);% save beta, SE

pRandomStat=NaN(size(randstat,1),2,numCompute);% Estimate,SEPred
GLMMmap.RandomEffects.Group=randstat.Group;
GLMMmap.RandomEffects.Level=randstat.Level;
GLMMmap.RandomEffects.Name=randstat.Name;
GLMMmap.RandomEffects.DF=randstat.DF;
GLMMmap.RandomEffects.DX=designMatrix(glmexample,'Random');
RandomStat=NaN(size(randstat,1),2,ySize2,xSize2);

FixMap2=FixMap(:,compuIdx);

if has_fsolve==1
    % parpool;
    try parpool(gridname);end
    parfor ip=1:numCompute
        fprintf('.')
        % construct lme table
        tbl = PredictorM;
        tbl.PixelIntensity = FixMap2(:,ip);
        % fit model
        glme = GeneralizedLinearMixedModel.fit(tbl,formula,varargin{:});
        % Ftest
        stats=anova(glme);
        % save result
        pAnovaFStat(:,ip)=double(stats(:,2));
        pCovb(:,:,ip)=glme.CoefficientCovariance;
        pDispersion(ip)=glme.Dispersion;
        % pMSEmatrix(ip)=glme.MSE;
        pSSEmatrix(ip)=glme.SSE;
        pSSTmatrix(ip)=glme.SST;
        pSSRmatrix(ip)=glme.SSR;
        pRsquaredMatrix(:,ip)=[glme.Rsquared.Ordinary glme.Rsquared.Adjusted];
        pModelCriterionMatrix(:,ip)=double(glme.ModelCriterion);% AIC,BIC,LogLikelihood,Deviance
        pCoefficientsMatrix(:,:,ip)= double(glme.Coefficients(:,[2 3]));% save beta, SE
        [~,~,stat]=randomEffects(glme);
        pRandomStat(:,:,ip)=double(stat(:,[4 5]));% Estimate,SEPred
    end
    % delete(gcp)
else
    h = waitbar(0,'Fitting GLMM, please wait...');
    for ip=1:numCompute
        waitbar(ip / numCompute)
        % construct lme table
        tbl = PredictorM;
        tbl.PixelIntensity = FixMap2(:,ip);
        % fit model
        glme = GeneralizedLinearMixedModel.fit(tbl,formula,varargin{:});
        % Ftest
        stats=anova(glme);
        % save result
        pAnovaFStat(:,ip)=double(stats(:,2));
        pCovb(:,:,ip)=glme.CoefficientCovariance;
        % pMSEmatrix(ip)=glme.MSE;
        pDispersion(ip)=glme.Dispersion;
        pSSEmatrix(ip)=glme.SSE;
        pSSTmatrix(ip)=glme.SST;
        pSSRmatrix(ip)=glme.SSR;
        pRsquaredMatrix(:,ip)=[glme.Rsquared.Ordinary glme.Rsquared.Adjusted];
        pModelCriterionMatrix(:,ip)=double(glme.ModelCriterion);% AIC,BIC,LogLikelihood,Deviance
        pCoefficientsMatrix(:,:,ip)= double(glme.Coefficients(:,[2 3]));% save beta, SE
        [~,~,stat]=randomEffects(glme);
        pRandomStat(:,:,ip)=double(stat(:,[4 5]));% Estimate,SEPred
    end
    close(h)
end

AnovaFStat(:,compuIdx)=pAnovaFStat;
Covb(:,:,compuIdx)=pCovb;
Dispersion(compuIdx)=pDispersion;
% MSEmatrix(compuIdx)=pMSEmatrix;
SSEmatrix(compuIdx)=pSSEmatrix;
SSTmatrix(compuIdx)=pSSTmatrix;
SSRmatrix(compuIdx)=pSSRmatrix;
RsquaredMatrix(:,compuIdx)=pRsquaredMatrix;
ModelCriterionMatrix(:,compuIdx)=pModelCriterionMatrix;
CoefficientsMatrix(:,:,compuIdx)=pCoefficientsMatrix;

GLMMmap.Anova.FStat=AnovaFStat;
GLMMmap.CoefficientCovariance=Covb;
% GLMMmap.MSE=MSEmatrix;
GLMMmap.Dispersion=Dispersion;
GLMMmap.SSE=SSEmatrix;
GLMMmap.SST=SSTmatrix;
GLMMmap.SSR=SSRmatrix;
GLMMmap.Rsquared=RsquaredMatrix;
GLMMmap.ModelCriterion=ModelCriterionMatrix;% AIC,BIC,LogLikelihood,Deviance
GLMMmap.Coefficients=CoefficientsMatrix;% save beta, SE

RandomStat(:,:,compuIdx)=pRandomStat;
GLMMmap.RandomEffects.RandomStat=RandomStat;
end