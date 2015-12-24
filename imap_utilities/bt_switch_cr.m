function [StatMap]=bt_switch_cr(StatMap,bootopt)
% bootstrap clustering test using other cluster criteria
% bootopt takes 1 cluster mass, 2 cluster size, 3 both, 4 cluster dense

mccopt=StatMap.mccopt;
statopt=StatMap.opt;
alpha=statopt.alpha;
c=statopt.c;
% h=statopt.h;
Pmap=StatMap.Pmap;
mapvalue=StatMap.map;

if strcmp(mccopt.methods,'bootstrap')==1
    ResampStat=StatMap.resampMat;
else
    error('ResampStat no found, please first perform bootstrap clustering test using imapLMMmcc.')
end

Pmask=zeros(size(Pmap));
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
    switch bootopt
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
StatMap.Pmask=logical(Pmask);