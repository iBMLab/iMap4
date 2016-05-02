function [RDM, stRDM, pMat, pstMat] = rdmfixmap(FixMap, Mask, CondiVec, SbjVec)
% compute representational dissimilarity matrix of smoothed fixation map
% basic on Mahalanobis distance.
DS         = cat(2,dummyvar(CondiVec),dummyvar(SbjVec));
[ds,ia,ic] = unique(DS,'rows');
unicd      = CondiVec(ia);

if size(ds,1)<size(DS,1) % compute the conditonal mean for each subject
    fixmap = zeros(size(ds,1), sum(Mask(:) == 1));
    for ist = 1:length(ia)
        idx           = ic==ist;
        fixmap(ist,:) = nanmean(FixMap(idx,Mask),1);
    end
else
    fixmap = FixMap(ia,Mask);
end

fixmap(isnan(fixmap))=0;


y     = pdist(fixmap,'mahalanobis');
stRDM = squareform(y);
end