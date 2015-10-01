function [CatePredictor,contrast]=CategoricalName(LMMmap)
% get catigorical predictor from model beta combination
%--------------------------------------------------------------------------
% Copyright (C) iMap Team 2015
tbl=LMMmap.Variables;
DesignMatrix=LMMmap.modelX;
if isa(tbl,'dataset')
    VarNames=tbl.Properties.VarNames;
elseif isa(tbl,'table')
    VarNames=tbl.Properties.VariableNames;
    tbl=table2dataset(tbl);
end
categypredi=LMMmap.VariableInfo.InModel&LMMmap.VariableInfo.IsCategorical;
categyprediname=VarNames(categypredi);
label=cell(length(tbl),sum(categypredi));
strlabel=cell(length(tbl),1);
for icc=1:sum(categypredi)
    label(:,icc)=cellstr(eval(['tbl.' categyprediname{icc}]));
end
for ii=1:length(label)
    strlabel{ii,:}=strjoin(label(ii,:),'_');
end
[CatePredictor,btmp]=unique(strlabel,'rows');
CateContrast=DesignMatrix(btmp,:);
CateContrast(CateContrast(:)~=0&CateContrast(:)~=1&CateContrast(:)~=-1)=0;
contrast=num2cell(CateContrast,2);
