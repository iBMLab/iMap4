function handles =  reset_variables(handles,step)

switch step
    
 case 1 % 'importdata' delete all

if isfield(handles,'data')
handles = rmfield(handles,{'data','colnames','colnames_original','data_original'});
end
if isfield(handles,'RawMap_estimated')
handles = rmfield(handles,{'RawMap_estimated','FixMap_estimated','Mask_estimated'});
end
if isfield(handles,'table')
set(handles.table, 'visible','off')
end

if isfield(handles,'table1')
    
set(handles.table1, 'visible','off')

end

if isfield(handles,'RawMap_estimated_scaled')
handles = rmfield(handles,{'RawMap_estimated_scaled','FixMap_estimated_scaled','Mask_estimated_scaled'});
end

if isfield(handles,'predictors_continuous')
handles = rmfield(handles,{'predictors_categorical','predictors_continuous','predictors_continuous_default','predictors_categorical_default'});
end
%remove trial method
if isfield(handles,'RawMap_single_trial_scaled')
handles = rmfield(handles,{'RawMap_single_trial_scaled','FixMap_single_trial_scaled','Mask_single_trial_scaled'});
end 


%remove scaled method
if isfield(handles,'FixMap_estimated_scaled')
handles = rmfield(handles,{'FixMap_estimated_scaled','Mask_estimated_scaled'});
end

case 2 %check_column
  handles.data = handles.data_original;
  handles.colnames = handles.colnames_original;
  
%remove trial method
if isfield(handles,'RawMap_single_trial_scaled')
handles = rmfield(handles,{'RawMap_single_trial_scaled','FixMap_single_trial_scaled','Mask_single_trial_scaled'});
end
  
  
if isfield(handles,'predictors_continuous')
handles = rmfield(handles,{'predictors_categorical','predictors_continuous','predictors_continuous_default','predictors_categorical_default'});
end  
case 3 % parameters
  
if isfield(handles,'predictors_continuous')
handles = rmfield(handles,{'predictors_categorical','predictors_continuous','predictors_continuous_default','predictors_categorical_default'});
end
%remove trial method
if isfield(handles,'RawMap_single_trial_scaled')
handles = rmfield(handles,{'RawMap_single_trial_scaled','FixMap_single_trial_scaled','Mask_single_trial_scaled'});
end

%remove estimated method
if isfield(handles,'RawMap_estimated')
handles = rmfield(handles,{'RawMap_estimated','FixMap_estimated','Mask_estimated'});
end

%remove estimated method
if isfield(handles,'RawMap_estimated_scaled')
handles = rmfield(handles,{'RawMap_estimated_scaled','FixMap_estimated_scaled','Mask_estimated_scaled'});
end

    
case 4 % predictors
if isfield(handles,'predictors_continuous')
handles = rmfield(handles,{'predictors_categorical','predictors_continuous','predictors_continuous_default','predictors_categorical_default'});
end
%remove trial method
if isfield(handles,'RawMap_single_trial_scaled')
handles = rmfield(handles,{'RawMap_single_trial_scaled','FixMap_single_trial_scaled','Mask_single_trial_scaled'});
end 


if isfield(handles,'FixMap_estimated')
handles = rmfield(handles,{'FixMap_estimated','Mask_estimated','RawMap_estimated'});
end

if isfield(handles,'FixMap_estimated_scaled')
handles = rmfield(handles,{'FixMap_estimated_scaled','Mask_estimated_scaled','RawMap_estimated_scaled'});
end 

       
end    
end