function PredictorM = nominal_predictorM(all_categorical_conditions,duration,PredictorM)

headers_predictorM = PredictorM.Properties.VarNames;

for h = 1:length(headers_predictorM)
    idx = strcmp(headers_predictorM(h),all_categorical_conditions)==1;
    if sum(idx)~=0
        PredictorM.(headers_predictorM{h}) = nominal(PredictorM.(headers_predictorM{h}));
    else
        if strcmp(headers_predictorM(h),duration)==0 && ~isnumeric(PredictorM.(headers_predictorM{h}))%duration already double
            PredictorM.(headers_predictorM{h}) = str2double(PredictorM.(headers_predictorM{h}));
        end
    end    
end
