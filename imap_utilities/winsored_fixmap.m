%--------------------------------------------------------------------------
% Copyright (C) iMap Team 2015
%% remove NaN
remv=find(isnan(FixMap(:,1,1)));
FixMap(remv,:,:)=[];
TblST(remv,:)=[];
%% robust option - winsorization on a subject level 
% the idea here is for each participant there is likely some outliner
% trials. This procedure is the most appropriate if you do single trial
% analysis.
sbjidx=unique(TblST.sbjST);
Ns=length(sbjidx);
fixmapSTw=zeros(size(FixMap));
figure;
for is2=1:Ns
    %%
    is=(is2);
    % figure;
    tmpsel=find(TblST.sbjST==sbjidx(is));
    tmpmap=FixMap(tmpsel,:,:);
%     subplot(1,2,1)
%     imagesc(squeeze(mean(tmpmap,1)));
%     axis equal off
    % winsoring
    tmpmaplow=prctile(tmpmap,5);
    tmpmaphigh=prctile(tmpmap,95);
    highmask=tmpmap>repmat(tmpmaphigh,size(tmpmap,1),1,1);
    lowmask=tmpmap<repmat(tmpmaplow,size(tmpmap,1),1,1);
    remainmask=(highmask+lowmask)==0;
    tmpmap1=bsxfun(@times,highmask,tmpmaphigh);
    tmpmap2=bsxfun(@times,lowmask,tmpmaplow);
    tmpmap3=tmpmap;
    tmpmap3(remainmask==0)=0;
    tmpmapnew=tmpmap1+tmpmap2+tmpmap3;
    fixmapSTw(tmpsel,:,:)=tmpmapnew;
%     subplot(1,2,2)
%     imagesc(squeeze(mean(tmpmapnew,1)));
%     axis equal off
    subplot(5,7,is2)
    imagesc(squeeze(mean(tmpmapnew,1)));
    axis equal off
    title(num2str(Tbl{is,2}))
end
