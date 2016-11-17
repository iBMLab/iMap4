function [contrast, hypothesis] = maincontrast(coefname,maplabel,tblvarname)
%--------------------------------------------------------------------------
% Copyright (C) iMap Team 2015

% first unify catigorical predictor in interactions
coefname2 = cell(size(coefname));
for icoef = 1:length(coefname)
    nametmp = coefname{icoef};
    textidx = strfind(nametmp,':');
    match1 = strcmp(nametmp,maplabel);
    if sum(match1) == 1
        coefname2{icoef} = nametmp;
    else
        if ~isempty(textidx)
            namepart = [];
            for itx = 1:(length(textidx)+1)
                if itx == 1
                    namepart2 = nametmp(1:(textidx(itx)-1));
                elseif itx == (length(textidx)+1)
                    namepart2 = nametmp((textidx(itx-1)+1):end);
                else
                    namepart2 = nametmp((textidx(itx-1)+1):(textidx(itx)-1));
                end
                while sum(strcmp(namepart2,tblvarname)) == 0
                    namepart2 = removelast(namepart2);
                end
                disp(namepart2)
                namepart  = [namepart ':' namepart2];
            end
            coefname2{icoef} = namepart(2:end);
        else
            if strcmp(nametmp,'(Intercept)')~=1
                while sum(strcmp(nametmp,tblvarname)) == 0
                    nametmp = removelast(nametmp);
                end
            end
            coefname2{icoef} = nametmp;
        end
    end
end

I = eye(length(coefname2));

contrast    = cell(length(maplabel),1);
% hypothesis  = cell(length(maplabel),1);
termkcols   = zeros(length(maplabel),length(coefname2));
for il = 1:length(maplabel)
    k                = strcmp(maplabel{il},coefname2);
    termkcols(il, k) = 1;
end
for il = 1:length(maplabel)
    contrast{il}     = I(logical(termkcols(il, :)), :);
end
hypothesis  = num2cell(zeros(length(contrast),1),2);
end

function a2 = removelast(textinput)
a2        = [];
[a1,b1]   = strtok(textinput,'_');
a2        = [a2 '_' a1];
while ~isempty(b1)
    [a1,b1] = strtok(b1(2:end),'_');
    if ~isempty(b1)
        a2 = [a2 '_' a1];
    end
end
a2        = a2(2:end);
end