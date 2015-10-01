function mccopt = select_option_randomfield(mccopt)

%-->if user select randomfield
%user input smoothing pixel size and scaling factor. and mccopt.sigma=smoothing*scale_factor (maybe this information should be save somewhere so the user dont need to memorized it... but for now let's just assume they have perfect memory)
dlg_title = 'sigma=smoothing*scaling factor';
prompt = {'smoothing pixel size:';'Scaling factor:          '};
formats(1,1).type   = 'edit';
formats(2,1).type   = 'edit';
defaultanswer = {'',''};
formats(1,1).size = 50;
formats(2,1).size = 50;
formats(1,1).format = 'integer';
formats(2,1).format = 'integer';

answer = inputsdlg(prompt, dlg_title, formats, defaultanswer);

if isempty (answer{1}) && isempty(answer{2})==0
smoothing = answer{1};
scale_factor = answer{2};
mccopt.sigma=smoothing*scale_factor;
else
   errordlg('One or both values are empty, unable to perform calculation') 
end

end