function modelling_part()

idx = modelling_analysis();
if idx ==1 % user did select a file
    %load  LMMmap
    [filename, pathname] = uigetfile('*.mat',' Select the FixMap & LMMmap from before to proceed','MultiSelect','on');
    if iscell(filename)==0 %user selects only multiple file
        errordlg('FixMap and LMMmap files need to be selected');
        uiwait(gcf)
    else
        fullpathname = strcat(pathname,filename); %user selects only one file
        if isempty(fullpathname)
            errordlg('No file selected')
        else
            
            for i = 1:length(filename)
                load(strcat(pathname,filename{i}))
            end
            
            if exist('LMMmap','var')&& exist('FixMap','var')
                model_plot(LMMmap,FixMap,pathname);
            else
                errordlg('The LMMmap and the FixMap are mandatory to proceed')
            end
        end
    end
end
end