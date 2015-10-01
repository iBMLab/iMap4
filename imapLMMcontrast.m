function [StatMap]=imapLMMcontrast(LMMmap,opt)
% Usage: [StatMap]=imapLMMcontrast(LMMmap,opt)
% output contrast and stat from original model without multiple comparison
% correction
%
% opt.type     - model/fixed/random/model beta/predictor beta
% opt.alpha    - default 0.05
% opt.c        - for coefficients and Catepredictors only, cell array
%               containing contrast vector/matrix
% opt.h        - for coefficients and Catepredictors only, cell array
%               containing hypothesis vector/matrix
% opt.onetail  - option to do onetail test, perform on two tail threshold
%               for convenience (alpha/2)
% opt.name     - for coefficients and Catepredictors only, name of each
%               contrast (for plotting)
%
% output structure format {opt} {Pmap} {Pmask} {F/Tmap} {betamap(optional)}
% {labels of the maps}
% see also imapLMM, imapLMMdisplay
%
% 2015-02-12 Junpeng Lao, University of Fribourg.
%
% Copyright (C) iMap Team 2015

StatMap=struct;
if ~isfield(opt,'alpha')
    alpha=0.05;
    opt.alpha=alpha;
else
    alpha=opt.alpha;
end
assert( isnumeric(alpha) & isreal(alpha) & isscalar(alpha) );
assert( alpha >= 0 & alpha <= 1 );

mask=isnan(LMMmap.MSE)==0;
nonnan=find(mask);

opt.type=lower(opt.type);
switch opt.type
    case 'model'% output model fitting and criterion map
        %%
        maplabel={'R2-Ordinary';'R2-Adjusted';'AIC';'BIC';'LogLikelihood';'Deviance'};
        mapvalue=[LMMmap.Rsquared;LMMmap.ModelCriterion];
        StatMap.label=maplabel;% name
        StatMap.map=mapvalue;% statvalue
    case 'fixed'% output Fvalue map, pvalue map, and mask according to MCC
        %%
        maplabel=LMMmap.Anova.Term;
        mapvalue=LMMmap.Anova.FStat;
        DF1=LMMmap.Anova.DF1;
        DF2=LMMmap.Anova.DF2;
        if strcmp(maplabel{1},'(Intercept)')==1
            maplabel(1)=[];
            mapvalue(1,:,:)=[];
            DF1(1)=[];
            DF2(1)=[];
        end
        % construct contrast cell (for resampling)
        coefname=LMMmap.CoefficientNames;
        [opt.c,opt.h]=maincontrast(coefname,maplabel);
        StatMap.label=maplabel;% name
        StatMap.map=mapvalue;% statvalue
        StatMap.df=[DF1 DF2];
        DF1tmp=(repmat(DF1,[1,length(mapvalue(:,:))]));
        DF2tmp=(repmat(DF2,[1,length(mapvalue(:,:))]));
        pmaptmp=1-fcdf(mapvalue(:,:),DF1tmp,DF2tmp);
        Pmap=reshape(pmaptmp,size(mapvalue));
        StatMap.Pmap=Pmap;
    case 'random'% output Fvalue map, pvalue map, and mask according to MCC
        %%
        betamap=squeeze(LMMmap.RandomEffects.RandomStat(:,1,:,:));
        SEtmp=squeeze(LMMmap.RandomEffects.RandomStat(:,2,:,:));
        df=LMMmap.RandomEffects.DF;
        DF=repmat(df,[1,size(betamap,2),size(betamap,3)]);
        maplabeltmp=strcat(cellstr(LMMmap.RandomEffects.Group), ...
            cellstr(LMMmap.RandomEffects.Level),'_', ...
            cellstr(LMMmap.RandomEffects.Name));
        StatMap.label=maplabeltmp;% name
        mapvalue=betamap./SEtmp;
        Pmap=1-tcdf(mapvalue,DF);
        StatMap.df=[ones(length(df),1) df];
        StatMap.map=mapvalue;% statvalue
        StatMap.beta=betamap;
        StatMap.Pmap=Pmap;
    case 'predictor beta'% output F/Tvalue map, beta map, pvalue map, and mask according to MCC
        %%
        if ~isfield(LMMmap,'SinglePred')
            error('You need to estimate the predictor beta in the imapLMM')
        else
            maplabeltmp=LMMmap.SinglePred.CatePredictor;
            betatmp=squeeze(LMMmap.SinglePred.beta(:,1,:,:));
            % SEtmp=squeeze(LMMmap.SinglePred.beta(:,2,:,:));
            covb=LMMmap.SinglePred.Covb;
            DF2tmp=LMMmap.modelDFE;
            if ~isfield(opt,'c')
                disp('Contrast not provided, output each predictor map')
                contrast=num2cell(diag(ones(size(maplabeltmp))),2);
                if ~isfield(opt,'h')
                    disp('Hypothesis not provided, contrast against 0')
                    hypothesis=num2cell(zeros(length(contrast),1),2);
                else
                    hypothesis=opt.h;
                end
                opt.name=maplabeltmp;
            else
                contrast=opt.c;
                if ~iscell(contrast)
                    contrast={contrast};
                end
                if ~isfield(opt,'h')
                    disp('Hypothesis not provided, contrast against 0')
                    hypothesis=num2cell(zeros(length(contrast),1),2);
                else
                    hypothesis=opt.h;
                end
                if ~iscell(hypothesis)
                    hypothesis={hypothesis};
                end
            end
            if length(contrast)~=length(hypothesis)
                error('contrast cell and hypothesis cell length mismatch')
            end
            if isfield(opt,'name')
                contrastlabel=opt.name;
            else
                tmpv=1:length(hypothesis);
                contrastlabel=strcat({'contrast'},cellstr(num2str(tmpv')));
            end
            % contrast
            DF2=repmat(DF2tmp,length(hypothesis),1);
            DF1=zeros(size(DF2));
            mapvalue=NaN(length(hypothesis),size(betatmp,2),size(betatmp,3));
            Pmap=mapvalue;
            betamap=mapvalue;
            betaCI=NaN(length(hypothesis),2,size(betatmp,2),size(betatmp,3));
            for ic=1:length(hypothesis)
                c=unique(contrast{ic},'rows');
                h=hypothesis{ic};
                DF1(ic)=rank(c);
                betaCItmp=NaN(2,size(betatmp,2),size(betatmp,3));
                covbtmp=covb(:,:,nonnan);
                proj=covproj(c,covbtmp); % (c*covb(:,:,ip)*c')
                cbeta_h=(c*betatmp(:,nonnan)-h)';
                if size(c,1)>1
                    betamap(ic,nonnan)=mean(cbeta_h,2);
                    invk1=zeros(size(proj));
                    for ip2=1:length(nonnan)
                        invk1(:,:,ip2)=inv(proj(:,:,ip2));
                    end
                else
                    betamap(ic,:)=c*betatmp(:,:);
                    delta = squeeze(tinv(1-alpha/2,DF2(ic)).*sqrt(proj))';
                    betaCItmp(:,nonnan)=[betamap(ic,nonnan)-delta;betamap(ic,nonnan)+delta];
                    invk1=proj.^-1;
                end
                
                Ftmp1=squeeze(full(covproj(cbeta_h,invk1)./DF1(ic)));
                ptmp1=1-fcdf(Ftmp1,DF1(ic).*ones(length(Ftmp1),1),DF2(ic).*ones(length(Ftmp1),1));
                if isfield(opt,'onetail')
                    if strcmp(opt.onetail,'>')==1
                        onetailtmp=logical(nansum(cbeta_h<0,2));
                        Ftmp1(onetailtmp)=0;
                        ptmp1(onetailtmp)=1;
                    elseif strcmp(opt.onetail,'<')==1
                        onetailtmp=logical(nansum(cbeta_h>0,2));
                        Ftmp1(onetailtmp)=0;
                        ptmp1(onetailtmp)=1;
                    end
                end
                mapvalue(ic,nonnan)=Ftmp1;
                Pmap(ic,nonnan)=ptmp1;
                betaCI(ic,:,:,:)=betaCItmp;
            end
            opt.c=contrast;
            opt.h=hypothesis;
            StatMap.label=contrastlabel;% name
            StatMap.map=mapvalue;% statvalue
            StatMap.beta=betamap;
            StatMap.betaCI=betaCI;
            StatMap.df=[DF1 DF2];
            StatMap.Pmap=Pmap;
        end
    case 'model beta'% output F/Tvalue map, beta map, pvalue map, and mask according to MCC
        %%
        maplabeltmp=LMMmap.CoefficientNames;
        betatmp=squeeze(LMMmap.Coefficients(:,1,:,:));
        % SEtmp=squeeze(LMMmap.Coefficients(:,2,:,:));
        covb=LMMmap.CoefficientCovariance;
        DF2tmp=LMMmap.modelDFE;
        if ~isfield(opt,'c')
            disp('Contrast not provided, output each predictor map')
            contrast=num2cell(diag(ones(size(maplabeltmp))),2);
            if ~isfield(opt,'h')
                disp('Hypothesis not provided, contrast against 0')
                hypothesis=num2cell(zeros(length(contrast),1),2);
            else
                hypothesis=opt.h;
            end
            opt.name=maplabeltmp;
        else
            contrast=opt.c;
            if ~iscell(contrast)
                contrast={contrast};
            end
            if ~isfield(opt,'h')
                disp('Hypothesis not provided, contrast against 0')
                hypothesis=num2cell(zeros(length(contrast),1),2);
            else
                hypothesis=opt.h;
            end
            if ~iscell(hypothesis)
                hypothesis={hypothesis};
            end
        end
        if length(contrast)~=length(hypothesis)
            error('contrast cell and hypothesis cell length mismatch')
        end
        if isfield(opt,'name')
            contrastlabel=opt.name;
        else
            tmpv=1:length(hypothesis);
            contrastlabel=strcat({'contrast'},cellstr(num2str(tmpv')));
        end
        % contrast
        DF2=repmat(DF2tmp,length(hypothesis),1);
        DF1=zeros(size(DF2));
        mapvalue=NaN(length(hypothesis),size(betatmp,2),size(betatmp,3));
        Pmap=mapvalue;
        betamap=mapvalue;
        betaCI=NaN(length(hypothesis),2,size(betatmp,2),size(betatmp,3));
        for ic=1:length(hypothesis)
            c=unique(contrast{ic},'rows');
            h=hypothesis{ic};
            DF1(ic)=rank(c);
            betaCItmp=NaN(2,size(betatmp,2),size(betatmp,3));
            covbtmp=covb(:,:,nonnan);
            proj=covproj(c,covbtmp); % (c*covb(:,:,ip)*c')            
            cbeta_h=(c*betatmp(:,nonnan)-h)';
            if size(c,1)>1
                betamap(ic,nonnan)=mean(cbeta_h,2);
                invk1=zeros(size(proj));
                for ip2=1:length(nonnan)
                    invk1(:,:,ip2)=inv(proj(:,:,ip2));
                end
            else
                betamap(ic,:)=c*betatmp(:,:);
                delta = squeeze(tinv(1-alpha/2,DF2(ic)).*sqrt(proj))';
                betaCItmp(:,nonnan)=[betamap(ic,nonnan)-delta;betamap(ic,nonnan)+delta];
                invk1=proj.^-1;
            end
            
            Ftmp1=squeeze(full(covproj(cbeta_h,invk1)./DF1(ic)));
            ptmp1=1-fcdf(Ftmp1,DF1(ic).*ones(length(Ftmp1),1),DF2(ic).*ones(length(Ftmp1),1));
            if isfield(opt,'onetail')
                if strcmp(opt.onetail,'>')==1
                    onetailtmp=logical(nansum(cbeta_h<0,2));
                    Ftmp1(onetailtmp)=0;
                    ptmp1(onetailtmp)=1;
                elseif strcmp(opt.onetail,'<')==1
                    onetailtmp=logical(nansum(cbeta_h>0,2));
                    Ftmp1(onetailtmp)=0;
                    ptmp1(onetailtmp)=1;
                end
            end
            mapvalue(ic,nonnan)=Ftmp1;
            Pmap(ic,nonnan)=ptmp1;
            betaCI(ic,:,:,:)=betaCItmp;
        end
        opt.c=contrast;
        opt.h=hypothesis;
        StatMap.label=contrastlabel;% name
        StatMap.map=mapvalue;% statvalue
        StatMap.beta=betamap;
        StatMap.betaCI=betaCI;
        StatMap.df=[DF1 DF2];
        StatMap.Pmap=Pmap;
    otherwise
        error('Unexpected test type. Please specify opt.type as one of the following: ''model'', ''fixed'', ''random'', ''predictor beta'', or ''model beta''.');
end
if strcmp(opt.type,'model')~=1
    StatMap.Pmask=Pmap<alpha;
end
StatMap.opt=opt;
end