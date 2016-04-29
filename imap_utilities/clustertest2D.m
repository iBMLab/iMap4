function Pmask = clustertest2D(Fvalue,Pmap,alpha,thresMass,thresSize,thresDens)
%--------------------------------------------------------------------------
% Copyright (C) iMap Team 2015

tmpMsk             = Pmap<alpha;
[L,num]            = bwlabel(tmpMsk,4);
if ~isempty(thresDens)
    for il = 1:num
        if nansum(nansum(Fvalue(L==il)))/sum(sum(L==il)) < thresDens
            tmpMsk(L==il) = 0;
        end
    end
else
    if ~isempty(thresMass) && ~isempty(thresSize)% using both cluster size and stat sum value as criterion
        for il = 1:num
            if sum(sum(L==il)) < thresSize || nansum(nansum(Fvalue(L==il))) < thresMass
                tmpMsk(L==il) = 0;
            end
        end
    elseif isempty(thresMass) && ~isempty(thresSize) % using only cluster size
        for il=1:num
            if sum(sum(L==il)) < thresSize
                tmpMsk(L==il) = 0;
            end
        end
    elseif ~isempty(thresMass) && isempty(thresSize) % using only cluster stat nansum value
        for il=1:num
            if nansum(nansum(Fvalue(L==il))) < thresMass
                tmpMsk(L==il) = 0;
            end
        end
    end
end
Pmask = tmpMsk > 0;
end