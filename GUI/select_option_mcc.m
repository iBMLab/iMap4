function mccopt = select_option_mcc()
mccopt = struct;
prompt = {'Method:';'tfce:     '};
name = 'MCC options';
formats = struct('type', {}, 'style', {}, 'items', {},'format', {}, 'limits', {}, 'size', {});
formats(1,1).type   = 'list';
formats(1,1).format = 'integer';
formats(1,1).size = 80;
formats(1,1).items  = {'fdr','bonferroni','randomfield','cluster','bootstrap','permutation'};
formats(1,1).style  = 'popupmenu';

formats(2,1).type   = 'list';
formats(2,1).size = 50;
formats(2,1).style  = 'popupmenu';
formats(2,1).items  = {'0','1'};

defaultanswer = {5,1};
%formats(2,1).size = 70;
answer = inputsdlg(prompt, name, formats, defaultanswer);
mccopt.methods = formats(1,1).items{answer{1}};
mccopt.tfce   = str2double(formats(2,1).items{answer{2}});
