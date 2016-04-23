function mccopt = select_option_permutation(mccopt)
% mccopt.bootgroup- grouping variable for bootstrap and permutation (to 
%                   keep group variance constant). Input must be a cell 
%                   specifying the grouping variables in the PredictorM
% mccopt.nboot    - number of resampling for bootstrap or permutation, default 1000

prompt = {'bootgroup:     ';'nboot:            '};
name = 'Permutation options';
formats = struct('type', {}, 'style', {}, 'items', {},'format', {}, 'limits', {}, 'size', {});

formats(1,1).type   = 'edit';
formats(1,1).size = 70;
%formats(2,1).format = 'integer';

formats(2,1).type   = 'edit';
formats(2,1).size = 70;
formats(2,1).format = 'integer';


defaultanswer = {'',1000};
answer = inputsdlg(prompt, name, formats, defaultanswer);
tmp   = answer{1}; %check format
if ~iscell(tmp);tmp2{1}=tmp;else tmp2=tmp; end;
mccopt.bootgroup= tmp2;
mccopt.nboot       = answer{2};

