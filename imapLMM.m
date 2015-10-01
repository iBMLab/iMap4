function [LMMmap,lmexample]=imapLMM(FixMap,PredictorM,Mask,opt,formula,varargin)
% Usage: [LMMmap,lmexample]=imapLMM(FixMap,PredictorM,Mask,opt,formula,varargin)
% Input:    FixMap - total number of trials * xSize * ySize
%       PredictorM - dataset format of condition Matrix, total number of
%                    trials * number of predictor. Categorical column must
%                    set to nominal
%             Mask - 2D mask, for reducing the number of computation
%              opt - structure. option to define parallel grid (opt.parallelname)
%                    and option to compute each single categorical condition
%                    beta (opt.singlepredi)
% formula/varargin - same as fitlme, type '>>help fitlme' for more
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
%
% Copyright (C) iMap Team 2015
%% check Matlab version and toolbox before start
v=ver;
if exist('LinearMixedModel')==0
    error('some crucial toolboxes are missing in your current version of Matlab for imapLMM to run')
end
has_fsolve=any(strcmp({v.Name},'Parallel Computing Toolbox'));
if has_fsolve==1
    if isfield(opt,'parallelname')==0
        gridname='local';
    else
        gridname=opt.parallelname;
    end
    disp('imapLMM will be computed in parallel mode, you machine might become quite slow during the computation')
end
if isfield(opt,'singlepredi')==1;
    singlep=opt.singlepredi;
    if singlep~=0
        disp('Single predictor beta for each categorical condition will be computed')
    end
end
% create save structure
LMMmap=struct;
LMMmap.runopt=opt;
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
    error('Please reform your input predictor matrix as dataset')
end
lmexample = LinearMixedModel.fit(tbl,formula,varargin{:});
[~,~,randstat]=randomEffects(lmexample);
time=toc;
disp(lmexample)
estimatime=time*numCompute*(2-(singlep==0))/60;
disp(['The total computation time is around ' num2str(estimatime) ... 
    ' minutes without parallel computation'])
% save information
LMMmap.VariableInfo=lmexample.VariableInfo;
LMMmap.Variables=tbl;
LMMmap.FitMethod=lmexample.FitMethod;
LMMmap.Formula=lmexample.Formula;
LMMmap.modelX=designMatrix(lmexample);
LMMmap.FitOptions=varargin;
LMMmap.modelDFE=lmexample.DFE;
LMMmap.CoefficientNames=lmexample.CoefficientNames;
% save Anova information
statexample=anova(lmexample);% no need to save pValue=1-fcdf(FStat,DF1,DF2);
LMMmap.Anova.Term=statexample.Term;
LMMmap.Anova.DF1=statexample.DF1;
LMMmap.Anova.DF2=statexample.DF2;

AnovaFStat=NaN(size(statexample.Term,1),ySize2,xSize2);
Covb=NaN([size(lmexample.CoefficientCovariance),ySize2,xSize2]);
MSEmatrix=NaN(ySize2,xSize2);
SSEmatrix=MSEmatrix;
SSTmatrix=MSEmatrix;
SSRmatrix=MSEmatrix;
RsquaredMatrix=NaN(2,ySize2,xSize2);
ModelCriterionMatrix=NaN(4,ySize2,xSize2);% AIC,BIC,LogLikelihood,Deviance
CoefficientsMatrix=NaN(length(lmexample.CoefficientNames),2,ySize2,xSize2);% save beta, SE
% for parallel computing
pAnovaFStat=NaN(size(statexample.Term,1),numCompute);
pCovb=NaN([size(lmexample.CoefficientCovariance),numCompute]);
pMSEmatrix=NaN(1,numCompute);
pSSEmatrix=pMSEmatrix;
pSSTmatrix=pMSEmatrix;
pSSRmatrix=pMSEmatrix;
pRsquaredMatrix=NaN(2,numCompute);
pModelCriterionMatrix=NaN(4,numCompute);% AIC,BIC,LogLikelihood,Deviance
pCoefficientsMatrix=NaN(length(lmexample.CoefficientNames),2,numCompute);% save beta, SE

if singlep~=0
    coefname=lmexample.CoefficientNames;
    categypredi=lmexample.VariableInfo.InModel&lmexample.VariableInfo.IsCategorical;
    categyprediname=VarNames(categypredi);
    % remove subject/grouping/random predictor if there is any
    exclu1=zeros(length(categyprediname),1);
    for icate=1:length(categyprediname)
        if isempty(strmatch(categyprediname{icate},coefname))
            exclu1(icate)=1;
        end
    end
    categyprediname(exclu1==1)=[];
    Ncatepred=length(categyprediname);
    if Ncatepred==0
        singlep=0;
    else
        label=cell(length(tbl),Ncatepred);
        strlabel=cell(length(tbl),1);
        for icc=1:Ncatepred
            label(:,icc)=cellstr(eval(['tbl.' categyprediname{icc}]));
        end
        for ii=1:size(label,1)
            strlabel{ii,:}=strjoin(label(ii,:),'_');
        end
        CatePredictor=unique(strlabel);CatePredictor=strcat({'c'},CatePredictor);
        ranprediname=cell(1,1);
        ranprediname{1}={'Intercept'};
        rangroupname{1}=lmexample.Formula.GroupingVariableNames{:};
        sbjvect=eval(['tbl.' char(rangroupname{:})]);
        DS=dummyvar(nominal(strlabel));
        CoefficientsMatrix2=NaN(length(CatePredictor),2,ySize2,xSize2);% save beta, SE
        Covb2=NaN([length(CatePredictor),length(CatePredictor),ySize2,xSize2]);
        pCoefficientsMatrix2=NaN(length(CatePredictor),2,numCompute);% save beta, SE
        pCovb2=NaN([length(CatePredictor),length(CatePredictor),numCompute]);
        LMMmap.SinglePred.CatePredictor=CatePredictor;
        LMMmap.SinglePred.DesignMatrix=DS;
    end
end
pRandomStat=NaN(size(randstat,1),2,numCompute);% Estimate,SEPred
LMMmap.RandomEffects.Group=randstat.Group;
LMMmap.RandomEffects.Level=randstat.Level;
LMMmap.RandomEffects.Name=randstat.Name;
LMMmap.RandomEffects.DF=randstat.DF;
LMMmap.RandomEffects.DX=designMatrix(lmexample,'Random');
RandomStat=NaN(size(randstat,1),2,ySize2,xSize2);

FixMap2=FixMap(:,compuIdx);

if has_fsolve==1
    % parpool;
    try parpool(gridname);end
    switch singlep
        case 0
            parfor ip=1:numCompute
                fprintf('.')
                % construct lme table
                tbl = PredictorM;
                tbl.PixelIntensity = FixMap2(:,ip);
                % fit model
                lme = LinearMixedModel.fit(tbl,formula,varargin{:});
                % Ftest
                stats=anova(lme);
                % save result
                pAnovaFStat(:,ip)=double(stats(:,2));
                pCovb(:,:,ip)=lme.CoefficientCovariance;
                pMSEmatrix(ip)=lme.MSE;
                pSSEmatrix(ip)=lme.SSE;
                pSSTmatrix(ip)=lme.SST;
                pSSRmatrix(ip)=lme.SSR;
                pRsquaredMatrix(:,ip)=[lme.Rsquared.Ordinary lme.Rsquared.Adjusted];
                pModelCriterionMatrix(:,ip)=double(lme.ModelCriterion);% AIC,BIC,LogLikelihood,Deviance
                pCoefficientsMatrix(:,:,ip)= double(lme.Coefficients(:,[2 3]));% save beta, SE
                [~,~,stat]=randomEffects(lme);
                pRandomStat(:,:,ip)=double(stat(:,[4 5]));% Estimate,SEPred
            end
        otherwise
            parfor ip=1:numCompute
                fprintf('.')
                % construct lme table
                tbl = PredictorM;
                tbl.PixelIntensity = FixMap2(:,ip);
                % fit model
                lme = LinearMixedModel.fit(tbl,formula,varargin{:});
                % Ftest
                stats=anova(lme);
                % save result
                pAnovaFStat(:,ip)=double(stats(:,2));
                pCovb(:,:,ip)=lme.CoefficientCovariance;
                pMSEmatrix(ip)=lme.MSE;
                pSSEmatrix(ip)=lme.SSE;
                pSSTmatrix(ip)=lme.SST;
                pSSRmatrix(ip)=lme.SSR;
                pRsquaredMatrix(:,ip)=[lme.Rsquared.Ordinary lme.Rsquared.Adjusted];
                pModelCriterionMatrix(:,ip)=double(lme.ModelCriterion);% AIC,BIC,LogLikelihood,Deviance
                pCoefficientsMatrix(:,:,ip)= double(lme.Coefficients(:,[2 3]));% save beta, SE
                lme2 = LinearMixedModel.fitmatrix(DS,tbl.PixelIntensity,ones(size(sbjvect)),sbjvect, ...
                    'FixedEffectPredictors',CatePredictor,'ResponseVarName','PixelIntensity','RandomEffectGroups',rangroupname{:}, ...
                    'RandomEffectPredictors',ranprediname);
                pCovb2(:,:,ip)=lme2.CoefficientCovariance;
                pCoefficientsMatrix2(:,:,ip)= double(lme2.Coefficients(:,[2 3]));% save beta, SE
                [~,~,stat]=randomEffects(lme);
                pRandomStat(:,:,ip)=double(stat(:,[4 5]));% Estimate,SEPred
            end
    end
    % delete(gcp)
else
    h = waitbar(0,'Fitting LMM, please wait...');
    for ip=1:numCompute
        waitbar(ip / numCompute)
        % construct lme table
        tbl = PredictorM;
        tbl.PixelIntensity = FixMap2(:,ip);
        % fit model
        lme = LinearMixedModel.fit(tbl,formula,varargin{:});
        % Ftest
        stats=anova(lme);
        % save result
        pAnovaFStat(:,ip)=double(stats(:,2));
        pCovb(:,:,ip)=lme.CoefficientCovariance;
        pMSEmatrix(ip)=lme.MSE;
        pSSEmatrix(ip)=lme.SSE;
        pSSTmatrix(ip)=lme.SST;
        pSSRmatrix(ip)=lme.SSR;
        pRsquaredMatrix(:,ip)=[lme.Rsquared.Ordinary lme.Rsquared.Adjusted];
        pModelCriterionMatrix(:,ip)=double(lme.ModelCriterion);% AIC,BIC,LogLikelihood,Deviance
        pCoefficientsMatrix(:,:,ip)= double(lme.Coefficients(:,[2 3]));% save beta, SE
        if singlep~=0
            lme2 = LinearMixedModel.fitmatrix(DS,tbl.PixelIntensity,ones(size(sbjvect)),sbjvect, ...
                'FixedEffectPredictors',CatePredictor,'ResponseVarName','PixelIntensity','RandomEffectGroups',rangroupname{:}, ...
                'RandomEffectPredictors',ranprediname);
            pCovb2(:,:,ip)=lme2.CoefficientCovariance;
            pCoefficientsMatrix2(:,:,ip)= double(lme2.Coefficients(:,[2 3]));% save beta, SE
        end
        [~,~,stat]=randomEffects(lme);
        pRandomStat(:,:,ip)=double(stat(:,[4 5]));% Estimate,SEPred
    end
    close(h)
end

AnovaFStat(:,compuIdx)=pAnovaFStat;
Covb(:,:,compuIdx)=pCovb;
MSEmatrix(compuIdx)=pMSEmatrix;
SSEmatrix(compuIdx)=pSSEmatrix;
SSTmatrix(compuIdx)=pSSTmatrix;
SSRmatrix(compuIdx)=pSSRmatrix;
RsquaredMatrix(:,compuIdx)=pRsquaredMatrix;
ModelCriterionMatrix(:,compuIdx)=pModelCriterionMatrix;
CoefficientsMatrix(:,:,compuIdx)=pCoefficientsMatrix;

LMMmap.Anova.FStat=AnovaFStat;
LMMmap.CoefficientCovariance=Covb;
LMMmap.MSE=MSEmatrix;
LMMmap.SSE=SSEmatrix;
LMMmap.SST=SSTmatrix;
LMMmap.SSR=SSRmatrix;
LMMmap.Rsquared=RsquaredMatrix;
LMMmap.ModelCriterion=ModelCriterionMatrix;% AIC,BIC,LogLikelihood,Deviance
LMMmap.Coefficients=CoefficientsMatrix;% save beta, SE
if singlep~=0
    Covb2(:,:,compuIdx)=pCovb2;
    CoefficientsMatrix2(:,:,compuIdx)= pCoefficientsMatrix2;
    
    LMMmap.SinglePred.beta=CoefficientsMatrix2;
    LMMmap.SinglePred.Covb=Covb2;
end
RandomStat(:,:,compuIdx)=pRandomStat;
LMMmap.RandomEffects.RandomStat=RandomStat;
end