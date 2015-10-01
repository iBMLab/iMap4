function iindex = retrieve_table_conditions(combMat,sTable,IDX_condition)

index =[];
if iscell(sTable(:,IDX_condition(1)))==1
    for i=1:size(combMat,2)
        index= [index,strcmp(sTable(:,IDX_condition(i)),combMat(1,i))];
    end
else
    for i=1:size(combMat,2)
        index= [index,sTable(:,IDX_condition(i))==combMat(1,i)];
    end
end

iindex = zeros(size(sTable,1),1);
idx = all(index,2)==1;
iindex(idx) = 1;
