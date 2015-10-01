function mccopt = select_option_bootstrap(mccopt)
% mccopt.bootopt  - 1 cluster mass, 2 cluster size, 3 both
% mccopt.bootgroup- grouping variable for bootstrap (to keep group variance
%                   constant). Input must be a cell specifying the grouping
%                   variables in the PredictorM
% mccopt.nboot    - number of resampling for bootstrap or permutation, default 1000

prompt = {'Bootopt:';'bootgroup:     ';'nboot:            '};
name = 'Bootstrap options';
formats = struct('type', {}, 'style', {}, 'items', {},'format', {}, 'limits', {}, 'size', {});
formats(1,1).type   = 'list';
formats(1,1).size = 100;
formats(1,1).items  = {'Cluster mass','Cluster size','Both','Dense'};
formats(1,1).style  = 'popupmenu';

formats(2,1).type   = 'edit';
formats(2,1).size = 70;
%formats(2,1).format = 'integer';

formats(3,1).type   = 'edit';
formats(3,1).size = 70;
formats(3,1).format = 'integer';


defaultanswer = {1,'',1000};
answer = inputsdlg(prompt, name, formats, defaultanswer);
mccopt.bootopt     = answer{1};
tmp   = answer{2}; %check format
if ~iscell(tmp);tmp2{1}=tmp;else tmp2=tmp; end;
mccopt.bootgroup= tmp2;
mccopt.nboot       = answer{3};

