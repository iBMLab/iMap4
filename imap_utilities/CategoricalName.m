function [CatePredictor,contrast,CateContrast] = CategoricalName(LMMmap)
% get catigorical predictor from model beta combination
%--------------------------------------------------------------------------
% Copyright (C) iMap Team 2015
tbl          = LMMmap.Variables;
DesignMatrix = LMMmap.modelX;
if     isa(tbl, 'dataset')
    VarNames = tbl.Properties.VarNames;
elseif isa(tbl, 'table')
    VarNames = tbl.Properties.VariableNames;
    tbl      = table2dataset(tbl);
end
categypredi     = LMMmap.VariableInfo.InModel ...
                & LMMmap.VariableInfo.IsCategorical;
categyprediname = VarNames(categypredi);
repredi         = LMMmap.Formula.GroupingVariableNames;
excld = zeros(length(categyprediname),length(repredi));
for irep = 1:length(repredi)
    excld(:,irep) = strcmp(categyprediname,repredi{irep});
end
categyprediname = categyprediname(sum(excld,2)==0);
label           = cell(length(tbl),length(categyprediname));
strlabel        = cell(length(tbl),1);

for icc = 1:length(categyprediname)
    label(:,icc)     = cellstr(eval(['tbl.' categyprediname{icc}]));
end
for ii = 1:length(label)
    strlabel{ii,:}   = strjoin(label(ii,:),'_');
end
[CatePredictor,btmp] = unique(strlabel,'rows');
CateContrast         = DesignMatrix(btmp,:);

idx                  = CateContrast(:) ~= 0 ...
                     & CateContrast(:) ~= 1 ...
                     & CateContrast(:) ~= -1;
CateContrast(idx)    = 0;

contrast             = num2cell(CateContrast,2);
