function [mask,cthred]=mcc2d(Statmap,pmap,alpha,method)
% multiple comparison correction base on the size of the clusters
% 
% naive clustering method
% input Statmap: F value or T value^2, 2 or 3 dimension data,
%                xSize*ySize*number of effects
%          pmap: p values
%         alpha: default .05
%        method: estimate all effect separetely -1; together -2

if nargin < 2
    error(message('Too Few Inputs'));
end
if nargin < 3
    alpha=.05;
    method=2;
end
if nargin < 4
    method=2;
end
% correction 
alpha=alpha/5;
correction=2.5;% 3
%%
if method==1
    mask=zeros(size(pmap));
    for iefft=1:size(Statmap,3)
        tempMAP=squeeze(Statmap(:,:,iefft));
        masktmp=zeros(size(pmap,1),size(pmap,2));
        pmask=squeeze(pmap(:,:,iefft))<alpha;
        [L, num] = bwlabel(pmask, 4);
        psum=zeros(1,num);
        for in=2:num
            psum(in)=sum(sum(L==in));
        end
        cthred(iefft)=max(psum)-(max(psum)-min(psum))*(.05*correction);
        indx=find(psum>cthred(iefft));
        for ind=1:length(indx)
            masktmp(L==indx(ind))=1;
        end
        mask(:,:,iefft)=logical(masktmp);
    end
else
    ik=0;
    mask=zeros(size(pmap));
    pmask=zeros(size(pmap));
    for iefft=1:size(Statmap,3)
        tempMAP=squeeze(Statmap(:,:,iefft));
        pmask(:,:,iefft)=squeeze(pmap(:,:,iefft))<alpha;
        [L, num] = bwlabel(squeeze(pmask(:,:,iefft)), 4);
        L2(:,:,iefft)=L;
        for in=2:num
            ik=ik+1;
            psum2(ik)=sum(sum(L==in));
        end
    end
    sortpvet=sort(psum2);
    cthred=sortpvet(ceil(length(psum2)*(1-.05*correction)));
    for iefft=1:size(Statmap,3)
        L=squeeze(L2(:,:,iefft));
        tempMAP=squeeze(Statmap(:,:,iefft));
        [L, num] = bwlabel(squeeze(pmask(:,:,iefft)), 4);
        psum=zeros(1,num);
        masktmp=zeros(size(pmap,1),size(pmap,2));
        for in=1:num
            psum(in)=sum(sum(L==in));
        end
        indx=find(psum>cthred);
        for ind=1:length(indx)
            masktmp(L==indx(ind))=1;
        end
        mask(:,:,iefft)=logical(masktmp);
    end
end