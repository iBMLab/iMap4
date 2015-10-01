function mccopt = select_option_fdr(mccopt)

%next window show them to chose 0 nonparametic or 1 parametic, result save in mccopt.parametic

prompt = {'Parametic:'};
dlg_title = 'Parametic Input';
formats = struct('type', {}, 'style', {}, 'items', {},'format', {}, 'limits', {}, 'size', {});
formats(1,1).type   = 'list';
formats(1,1).style  = 'popupmenu';
formats(1,1).items  = {'Non-parametic','Parametic'};
defaultanswer = {1};
formats(1,1).size = 100;
answer = inputsdlg(prompt, dlg_title, formats, defaultanswer);

if answer{1} == 1 % user choose first option -> nonparametic
mccopt.parametic = 0;
else
mccopt.parametic = 1; % user choose parametic

end
