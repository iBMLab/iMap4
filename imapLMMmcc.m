function [StatMap_c]=imapLMMmcc(StatMap,LMMmap,mccopt,varargin)
% {MCC related field}
% mccopt.methods  - fdr/bonferroni/randomfield/cluster/bootstrap/permutation
% mccopt.bootopt  - 1 cluster mass, 2 cluster size, 3 both cluster mass and
%                   size, 4 cluster dense 
% mccopt.bootgroup- grouping variable for bootstrap and permutation (to 
%                   keep group variance constant). Input must be a cell 
%                   specifying a Group variables in the PredictorM
% mccopt.sbjvec   - subject vector for bootstrap. Input must be a cell 
%                   specifying a Group variables in the PredictorM. This is 
%                   important when there are multiple grouping variables  
%                   in the mixed model such as (1|subject) + (1|stimuli)
% mccopt.nboot    - number of resampling for bootstrap or permutation
% mccopt.sigma    - smoothing parameter (for Random field test)
% mccopt.clustSize- cluster size threshold (for cluster test)
% mccopt.clustVal - cluster value threshold (for cluster test)
% mccopt.parametic- for FDR
% mccopt.tfce     - signal enhancement base on Threshold-free cluster
%                   enhancement developed by Smith & Nichols, 2009
% varargin        - replace it with [FixMap] for resampling algorithm.
% 
% 2015-02-12 Junpeng Lao, University of Fribourg.
%--------------------------------------------------------------------------
% Copyright (C) iMap Team 2015

warning('off','all');
if nargin > 3
    FixMap=varargin{1};
end
statopt        = StatMap.opt;
alpha          = statopt.alpha;
Pmap           = StatMap.Pmap;
mapvalue       = StatMap.map;
StatMap.Pmask0 = StatMap.Pmask;
% clear old fields
StatMap.Pmask  = [];
mask   = isnan(LMMmap.MSE)==0;
tbl    = LMMmap.Variables;
nonnan = find(mask);

%% multiple comparison correction
if strcmp(statopt.type,'model')~=1
    Pmask          = zeros(size(Pmap));
    mccopt.methods = lower(mccopt.methods);
    switch mccopt.methods
        case 'fdr'
            %%
            pid = zeros(1,size(Pmask,1));
            if isfield(mccopt,'parametic')==1
                for imask=1:size(Pmask,1)
                    tmpMsk           = squeeze(Pmap(imask,:,:));
                    pid(imask)       = fdr(tmpMsk(isnan(tmpMsk)==0), alpha,'parametic');
                    Pmask(imask,:,:) = tmpMsk<pid(imask);
                end
            else
                for imask=1:size(Pmask,1)
                    tmpMsk           = squeeze(Pmap(imask,:,:));
                    pid(imask)       = fdr(tmpMsk(isnan(tmpMsk)==0), alpha);
                    Pmask(imask,:,:) = tmpMsk<pid(imask);
                end
            end
            % tfec
            if isfield(mccopt,'tfce')==1
                if mccopt.tfce==1
                    mapvaluetmp      = StatMap.map;
                    mapvalue2        = tfce2d(permute(mapvaluetmp,[2,3,1]));
                    mapvalue         = permute(mapvalue2,[3,1,2]);
                    StatMap.map      = mapvalue;
                end
            end
            mccopt.fdrThreshold = pid;
        case 'bonferroni'
            %%
            alpha2                     = alpha/length(nonnan);
            Pmask                      = Pmap<alpha2;
            mccopt.BonferroniThreshold = alpha2;
            % tfec
            if isfield(mccopt,'tfce')==1
                if mccopt.tfce==1
                    mapvaluetmp = StatMap.map;
                    mapvalue2   = tfce2d(permute(mapvaluetmp,[2,3,1]));
                    mapvalue    = permute(mapvalue2,[3,1,2]);
                    StatMap.map = mapvalue;
                end
            end
        case 'randomfield'
            %%
            % for details see Chauvin, A., Worsley, K. J., Schyns, P. G., Arguin, M. &
            % Gosselin, F. (2004).  A sensitive statistical test for smooth
            % classification images.
            pixsearchspace = mask;
            % calculation of the significant values for single matrix
            sigma          = mccopt.sigma;
            FWHM           = sigma * sqrt(8*log(2));% computes the full width half maximum
            [volumes,N]    = CiVol(pixsearchspace,2); % (Worsley et al. 1996, HBM) computes the intrinsic volumes
            tP             = zeros(1,size(Pmask,1));
            for imask=1:size(Pmask,1)
                Fvalue     = squeeze(mapvalue(imask,:,:));
                if strcmp(statopt.type,'random');Fvalue = Fvalue.^2;end
                tP(imask)        = stat_threshold(volumes, N,FWHM,StatMap.df(imask,:),alpha);
                Pmask(imask,:,:) = Fvalue>tP(imask);
            end
            mccopt.randomfieldThreshold = tP;
            % tfec
            if isfield(mccopt,'tfce')==1
                if mccopt.tfce==1
                    mapvaluetmp = StatMap.map;
                    mapvalue2   = tfce2d(permute(mapvaluetmp,[2,3,1]));
                    mapvalue    = permute(mapvalue2,[3,1,2]);
                    StatMap.map = mapvalue;
                end
            end
        case 'cluster'% cluster size and stat value sum within cluster could be estimated via resampling (bootstrap or jackknife)
            %%
            if isfield(mccopt,'tfce')==1
                if mccopt.tfce==1
                    mapvaluetmp = StatMap.map;
                    mapvalue2   = tfce2d(permute(mapvaluetmp,[2,3,1]));
                    mapvalue    = permute(mapvalue2,[3,1,2]);
                    StatMap.map = mapvalue;
                end
            end
            if isfield(mccopt,'clustSize') && isfield(mccopt,'clustVal') % using both cluster size and stat sum value as criterion
                clustersize  = mccopt.clustSize;
                clustervalue = mccopt.clustVal;
                if length(clustersize)~=size(Pmap,1) || length(clustervalue)~=size(Pmap,1)
                    warning('cluster threshold missmatch, using the first element of the threshold vector')
                    clustersize2  = clustersize(1)*ones(1,size(Pmap,1));
                    clustervalue2 = clustervalue(1)*ones(1,size(Pmap,1));
                else
                    clustersize2  = clustersize;
                    clustervalue2 = clustervalue;
                end
                for imask=1:size(Pmask,1)
                    Fvalue = squeeze(mapvalue(imask,:,:));
                    pvalue = squeeze(Pmap(imask,:,:));
                    if strcmp(statopt.type,'random');Fvalue = Fvalue.^2;end
                    Pmask(imask,:,:) = clustertest2D(Fvalue,pvalue,alpha,clustervalue2(imask),clustersize2(imask),[]);
                end
            elseif isfield(mccopt,'clustSize') && ~isfield(mccopt,'clustVal') % using only cluster size
                clustersize = mccopt.clustSize;
                if length(clustersize)~=size(Pmap,1)
                    warning('cluster threshold missmatch, using the first element of the threshold vector')
                    clustersize2 = clustersize(1)*ones(1,size(Pmap,1));
                else
                    clustersize2 = clustersize;
                end
                for imask=1:size(Pmask,1)
                    Fvalue = squeeze(mapvalue(imask,:,:));
                    pvalue = squeeze(Pmap(imask,:,:));
                    if strcmp(statopt.type,'random');Fvalue = Fvalue.^2;end
                    Pmask(imask,:,:) = clustertest2D(Fvalue,pvalue,alpha,[],clustersize2(imask),[]);
                end
            elseif ~isfield(mccopt,'clustSize') && isfield(mccopt,'clustVal') % using only cluster stat sum value
                clustervalue = mccopt.clustVal;
                if length(clustervalue)~=size(Pmap,1)
                    warning('cluster threshold missmatch, using the first element of the threshold vector')
                    clustervalue2 = clustervalue(1)*ones(1,size(Pmap,1));
                else
                    clustervalue2 = clustervalue;
                end
                for imask=1:size(Pmask,1)
                    Fvalue = squeeze(mapvalue(imask,:,:));
                    pvalue = squeeze(Pmap(imask,:,:));
                    if strcmp(statopt.type,'random');Fvalue=Fvalue.^2;end
                    Pmask(imask,:,:) = clustertest2D(Fvalue,pvalue,alpha,clustervalue2(imask),[],[]);
                end
            end
        case 'bootstrap'
            %%
            if isfield(mccopt,'sbjvec')==1
                sbjvectmp = mccopt.sbjvec{1};
                sbjvec    = eval(['tbl.' sbjvectmp]);
            else
                sbjvec    = [];
            end
            grouping=nominal(ones(size(tbl,1),1));
            if isfield(mccopt,'bootgroup')==1
                Ng=length(mccopt.bootgroup);
                grouping1=nominal(zeros(size(tbl,1),Ng));
                for ig=1:Ng
                    groupingtmp=mccopt.bootgroup{ig};
                    if ~isempty(groupingtmp)
                    grouptmp=eval(['tbl.' groupingtmp]);
                    grouping1(:,ig)=nominal(grouptmp);
                    end
                end
                if Ng>1
                    for ii=1:length(grouping1)
                        tmp=char(grouping1(ii,:));
                        tmp=tmp';
                        grouping(ii,:)=cellstr(tmp(:)');
                    end
                else
                    grouping=grouping1;
                end
            end
            nboot=mccopt.nboot;
            c=statopt.c;
            h=statopt.h;
            [ResampStat]=imapLMMresample(FixMap,LMMmap,c,h,statopt.type,'bootstrap',nboot,grouping,1,sbjvec);
            % tfce on orignial
            if isfield(mccopt,'tfce')==1
                if mccopt.tfce==1
                    mapvaluetmp=StatMap.map;
                    mapvalue2=tfce2d(permute(mapvaluetmp,[2,3,1]));
                    mapvalue=permute(mapvalue2,[3,1,2]);
                    StatMap.map=mapvalue;
                end
            end
            bootclustdist=zeros(length(c),nboot,3);
            for ic=1:length(c)
                Fboot=squeeze(ResampStat.resFvalue(:,ic,:,:));
                Pboot=squeeze(ResampStat.resPvalue(:,ic,:,:));
                Betabt=squeeze(ResampStat.resBeta(:,ic,:,:));
                % one tail test
                if isfield(statopt,'onetail')
                    if strcmp(statopt.onetail,'>')==1
                        Fboot(Betabt<0)=0;
                        Pboot(Betabt<0)=1;
                    elseif strcmp(statopt.onetail,'<')==1
                        Fboot(Betabt>0)=0;
                        Pboot(Betabt>0)=1;
                    end
                end
                % tfce
                if isfield(mccopt,'tfce')==1
                    if mccopt.tfce==1
                        mapvalue2=tfce2d(permute(Fboot,[2,3,1]));
                        Fboot=permute(mapvalue2,[3,1,2]);
                    end
                end
                % output estimated bootstrap distribution (cluster test)
                maxclust=zeros(nboot,3);
                for ib=1:nboot
                    Ftmp=squeeze(Fboot(ib,:,:));
                    [bmasktmp,bnum]=bwlabel(squeeze(Pboot(ib,:,:))<alpha);
                    if bnum>0
                        maxtmp=zeros(bnum,3);
                        for icluster=1:bnum
                            maxtmp(icluster,1)=nansum(nansum(Ftmp(bmasktmp==icluster)));
                            maxtmp(icluster,2)=sum(sum(bmasktmp==icluster));
                            maxtmp(icluster,3)=maxtmp(icluster,1)./maxtmp(icluster,2);
                        end
                        maxclust(ib,1)=max(maxtmp(:,1));
                        maxclust(ib,2)=max(maxtmp(:,2));
                        maxclust(ib,3)=max(maxtmp(:,3));
                    end
                end
                distmp=sort(maxclust,1);
                bootclustdist(ic,:,:)=distmp;
                % output cluster mass and size threshold at alpha
                clthres=distmp(round(nboot*(1-alpha)),:);
                % output new Pmask: 1 cluster mass, 2 cluster size, 3 both,
                % 4 cluster dense
                switch mccopt.bootopt
                    case 1
                        Pmask(ic,:,:)=clustertest2D(squeeze(mapvalue(ic,:,:)), ...
                            squeeze(Pmap(ic,:,:)),alpha,clthres(1),[],[]);
                    case 2
                        Pmask(ic,:,:)=clustertest2D(squeeze(mapvalue(ic,:,:)), ...
                            squeeze(Pmap(ic,:,:)),alpha,[],clthres(2),[]);
                    case 3
                        Pmask(ic,:,:)=clustertest2D(squeeze(mapvalue(ic,:,:)), ...
                            squeeze(Pmap(ic,:,:)),alpha,clthres(1),clthres(2),[]);
                    case 4
                        Pmask(ic,:,:)=clustertest2D(squeeze(mapvalue(ic,:,:)), ...
                            squeeze(Pmap(ic,:,:)),alpha,[],[],clthres(3));
                end
            end
            StatMap.bootdist=bootclustdist;
            StatMap.resampMat=ResampStat;
        case 'permutation'
            %%
            grouping=nominal(ones(size(tbl,1),1));
            if isfield(mccopt,'bootgroup')==1
                Ng=length(mccopt.bootgroup);
                grouping1=nominal(zeros(size(tbl,1),Ng));
                for ig=1:Ng
                    groupingtmp=mccopt.bootgroup{ig};
                    if ~isempty(groupingtmp)
                    grouptmp=eval(['tbl.' groupingtmp]);
                    grouping1(:,ig)=nominal(grouptmp);
                    end
                end
                if Ng>1
                    for ii=1:length(grouping1)
                        tmp=char(grouping1(ii,:));
                        tmp=tmp';
                        grouping(ii,:)=cellstr(tmp(:)');
                    end
                else
                    grouping=grouping1;
                end
            end
            nboot=mccopt.nboot;
            c=statopt.c;
            h=statopt.h;
            [ResampStat]=imapLMMresample(FixMap,LMMmap,c,h,statopt.type,'permutation',nboot,grouping,1);
            % tfce on orignial
            if isfield(mccopt,'tfce')==1
                if mccopt.tfce==1
                    mapvaluetmp=StatMap.map;
                    mapvalue2=tfce2d(permute(mapvaluetmp,[2,3,1]));
                    mapvalue=permute(mapvalue2,[3,1,2]);
                    StatMap.map=mapvalue;
                    % mapParti=ResampStat.Forg;
                    % mapvalue2=tfce2d(permute(mapvaluetmp,[2,3,1]));
                    % mapvalue=permute(mapvalue2,[3,1,2]);
                    % ResampStat.Forg=mapParti;
                end
            end
            Pmapnew=NaN(size(Pmap));
            Pmapnew2=NaN(size(Pmap));
            StatMap.origPmap=Pmap;
            for ic=1:length(c)
                Fboot=squeeze(ResampStat.resFvalue(:,ic,:,:));
                
                % tfce
                if isfield(mccopt,'tfce')==1
                    if mccopt.tfce==1
                        mapvalue2=tfce2d(permute(Fboot,[2,3,1]));
                        Fboot=permute(mapvalue2,[3,1,2]);
                    end
                end
                
                % output new Pmask
                origFmap=squeeze(StatMap.map(ic,:,:));
                orFmat=permute(repmat(origFmap,[1,1,nboot]),[3,1,2]);
                Fboot2=repmat(max(Fboot(:,:),[],2),[1,size(Fboot,2),size(Fboot,3)]);
                Pmapnew(ic,:,:)=(sum(Fboot>=orFmat,1)+1)./(nboot+1);
                Pmapnew(ic,isnan(origFmap))=1;
                Pmapnew2(ic,:,:)=(sum(Fboot2>=orFmat,1)+1)./(nboot+1);
                Pmapnew2(ic,isnan(origFmap))=1;
                
                % one tail test
                if isfield(statopt,'onetail')
                    if strcmp(statopt.onetail,'>')==1
                        Pmapnew(StatMap.beta<0)=1;
                        Pmapnew2(StatMap.beta<0)=1;
                    elseif strcmp(statopt.onetail,'<')==1
                        Pmapnew(StatMap.beta>0)=1;
                        Pmapnew2(StatMap.beta>0)=1;
                    end
                end
                Pmasktmp=squeeze(Pmapnew2(ic,:,:)<alpha);
                Pmask(ic,:,:)=Pmasktmp;
            end
            StatMap.Pmap=Pmapnew;
            StatMap.Pmapfwer=Pmapnew2;
            StatMap.resampMat=ResampStat;
        otherwise
            error('Unexpected MC correction type. Please specify mccopt.method as one of the following: ''FDR'', ''Bonferroni'', ''Randomfield'', ''cluster'', ''bootstrap'' or ''permutation''.');
    end
    StatMap.Pmask=logical(Pmask);
    StatMap.mccopt=mccopt;
else
    warning('No multiple comparsion correction is available for model fitting criteria')
end

% save data
StatMap_c=StatMap;

warning('on','all');
end
