function [RDM, stRDM, unicd] = rdmfixmap(FixMap, Mask, CondiVec, SbjVec, varargin)
% compute representational dissimilarity matrix of smoothed fixation map
% basic on Mahalanobis distance.
% Disclaimer: the multivariate distance is computed on the PCA component
% and scaled by the % of variance explained. This approach is not yet
% peer-reviewed but seems to work empirically (in our simulation)
%--------------------------------------------------------------------------
% Copyright (C) iMap Team 2016

if nargin > 4
    plotopt = varargin{1};
else
    plotopt = 1;
end    
CondiVec   = categorical(CondiVec);
SbjVec     = categorical(SbjVec);
% normalized map
for i = 1:size(FixMap,1)
    tmpmap         = FixMap (i,:,:);
    tmpnorm        = tmpmap-mean(tmpmap(:))./std(tmpmap(:));
    FixMap (i,:,:) = tmpnorm;
end
DS         = cat(2,dummyvar(CondiVec),dummyvar(SbjVec));
[ds,ia,ic] = unique(DS,'rows');
cv2        = CondiVec(ia);
unicd      = categories(cv2);
Nc         = length(unicd);

if size(ds,1)<size(DS,1) % compute the conditonal mean for each subject
    fixmap = zeros(size(ds,1), sum(Mask(:) == 1));
    for ist = 1:length(ia)
        idx           = ic==ist;
        fixmap(ist,:) = nanmean(FixMap(idx,Mask),1);
    end
else
    fixmap = FixMap(ia,Mask);
end

fixmap(isnan(fixmap)) = 0;
% fixmap                = normalizeX(fixmap);
%%
[~,mapscore,~,~,explained,~] = pca(fixmap);


list      = cumsum(explained) <= 99;
weight    = diag(explained(list));

mapreduce = mapscore(:,list);
y         = pdist(mapreduce,'mahalanobis');
stRDM     = squareform(y);

mu        = mean(mapreduce,1);
sd        = nancov (mapreduce,1);
RDM       = NaN(Nc,Nc);
% L         = length(CondiVec);

for ic1 = 1:Nc
    map1         = mapreduce( cv2==unicd(ic1) , :);
    mu1          = mean(map1, 1);
    sd1          = nancov(map1 ,1);
    for ic2 = 1:Nc
        if ic1 == ic2
            mu2  = mu;
            % C    = sd1;
        else
            map2 = mapreduce( cv2==unicd(ic2) , :);
            mu2  = mean(map2, 1);
            sd2  = nancov(map2, 1);
            % C    = (sd1+sd2)/2;
            
            % dist = MIim (map1, map2, L);
        end
        C        = sd;
        % C        = sd*weight;
        
        % dist     = sqrt((mu1-mu2)*pinv(C)*(mu1-mu2)');
        % dist     = pdist2(mu1,mu2,'seuclidean',explained(list));
        % dist     = pdist2(mu1,mu2,'corr');
        
        % Multivariance distance (Mahalanobis)
        % M1 (uncomment the below one)
        % dist     = pdist2(mu1,mu2,'mahalanobis',C);

        % Multivariance distance scaled by variance explain, closely related to correlation
        % M2 (uncomment the below two)
        diffmu   = mu1-mu2;
        dist     = sqrt(sum((diffmu'.*(diag(C).^-1).*diffmu').*explained(list)));
        
        RDM (ic1,ic2) = dist;
    end
end
% RDM = RDM./max(RDM(:));
if plotopt
    % display output
    scrsz=get(0,'ScreenSize');% get screen size for output display
    figure('Numbertitle','off','Name',...
        'Representational Dissimilarity Matrix (value shows multivariate distance)',...
        'Position',[1 1 scrsz(3) scrsz(4)]);
    subplot(1,2,1)
    imagesc(stRDM);
    title('stRDM')
    axis square off;
    subplot(1,2,2)
    imsqrmat(RDM, unicd);
    title('RDM')
end
end

%---------------------------------------------
% Normalize the data matrices X and Y to have unit norm
% this part of the codes is from Matlab pdist2 function
function [X, flag] = normalizeX(X)
Xmax = max(abs(X),[],2);
X2 = bsxfun(@rdivide,X,Xmax);
Xnorm = sqrt(sum(X2.^2, 2));

% Find out points for which distance cannot be computed.

% The norm will be NaN for rows that are all zeros, fix that for the test
% below.
Xnorm(Xmax==0) = 0;

% The norm will be NaN for rows of X that have any +/-Inf. Those should be
% Inf, but leave them as is so those rows will not affect the test below.
% The points can't be normalized, so any distances from them will be NaN
% anyway.

% Find points that are effectively zero relative to the point with largest norm.
flag =  any(Xnorm <= eps(max(Xnorm)));
Xnorm = Xnorm .* Xmax;
X = bsxfun(@rdivide,X,Xnorm);
end

%---------------------------------------------
% Mutual Information
function I = MIim (A,B,L)
na    = hist(A(:),L);
na    = na/sum(na);
nb    = hist(B(:),L);
nb    = nb/sum(nb);
papb  = na'*nb;

ma=min(A(:));
MA=max(A(:));
mb=min(B(:));
MB=max(B(:));

% Scale and round to fit in {0,...,L-1}
A=round((A-ma)*(L-1)/(MA-ma+eps));
B=round((B-mb)*(L-1)/(MB-mb+eps));
n2=zeros(L);
x=0:L-1;
for i=0:L-1
    n2(i+1,:) = histc(B(A==i),x,1);
end
pab   = n2/sum(n2(:));

I     = find(papb(:)>1e-12 & pab(:)>1e-12); % function support
y     = pab(I).*log2(pab(I)./papb(I));

I     = sum(y);
end
