function StatMap_c=multiple_comparison_tab(LMMmap,FixMap,StatMap,pathname)

if strcmp(StatMap.opt.type,'model')==1
    errordlg('There is no multiple comparison correction on model fitting')
    StatMap_c=[];
    return
    uiwait(gcf)
end

mccopt = select_option_mcc();

switch mccopt.methods
    case 'bootstrap'
        mccopt = select_option_bootstrap(mccopt);
    case 'fdr'
        mccopt = select_option_fdr(mccopt);
    case 'bonferroni'
    case 'randomfield'
        mccopt = select_option_randomfield(mccopt);
    case 'cluster'
        mccopt = select_option_cluster(mccopt);
    case 'permutation'
        %mccopt.permute   = str2double(options{2}); %permutation test type. 1 for pixel wise p-value, 2 for cluster wise p-value
        errordlg('permutation still under construnction');
end

%mccopt.methods   = options{1};             %fdr/bonferroni/randomfield/cluster/bootstrap/permutation
%mccopt.tfce    = str2double(options{2}); %signal enhancement base on Threshold-free cluster
%mccopt.bootopt   = str2double(options{3}); %1 cluster mass, 2 cluster size, 3 both
%mccopt.bootgroup = str2double(options{4}); %grouping variable for bootstrap (to keep group variance  %constant). Input must be a cell specifying the grouping variables in the PredictorM
%mmcopt.nboot = str2double(options{5});     %number of resampling for bootstrap or permutation
%mccopt.sigma = str2double(options{6});     %smoothing parameter (for Random field test)
%mccopt.clustSize = str2double(options{7}); %cluster size threshold (for cluster test)
%mccopt.clustVal  = str2double(options{8}); %cluster value threshold (for cluster test)
%mccopt.parametic = str2double(options{9}); %for FDR

%compute StatMap_c,
StatMap_c=imapLMMmcc(StatMap,LMMmap,mccopt,FixMap);
prompt = {'Please rename the StatMap',};
dlg_title = 'Saving StatMap';
num_lines = 1;
if strcmp(StatMap_c.opt.type,'fixed')
def = {'StatMap_ANOVA_mcc_'};
elseif strcmp(StatMap_c.opt.type,'predictor beta')
def = {'StatMap_Contrast_mcc_'};
end
answer = [];
while isempty(answer)
    answer = inputdlg(prompt,dlg_title,num_lines,def);
end
if strcmp(mccopt.methods,'bootstrap')
    %then ask whether visualize map
    answer_user2 = questdlg('Would you like to remove the bootstrap distribution? (removing the resmapling subfield could speed up the saveing and loading time)','Remove bootstrap distribution','yes','no','no');
    if strcmp(answer_user2,'yes')
    % rm bootstrap field to save space and saving time.
    StatMap_c=rmfield(StatMap_c,'resampMat');
    end
end
save(strcat(pathname,answer{1}), 'StatMap_c','-v7.3');

%then ask whether visualize map
answer_user = questdlg('Would you like to visualize the map?','Visualisation of the map','yes','no','no');

if strcmp(answer_user,'yes')
    % output figure;
    %option for LMMdisplay
    %normalized,backgroundfile,colourmap,colormaprange,distplot
    prompt={'Normalize map values (1 to normalize all maps)','Background file','Colormap (as being used in Matlab)','map range for optimal display (0 to 1)','map Value distribution (1 to display)'};
    name='Options for Display';
    numlines = [1 50;1 50;1 50;1 50;1 50];
    defaultanswer={'','','''jet''','',''};
    options=inputdlg(prompt,name,numlines,defaultanswer);
    normalized      = str2double(options{1});
    cmap            = options{3};
    backgroundfile  = char(options{2});
    colormaprange   = str2double(options{4});
    distplot        = str2double(options{5});
    
    %output figure;
    % imapLMMdisplay(StatMap,0)
    imapLMMdisplay(StatMap_c,normalized,backgroundfile(2:end-1),cmap(2:end-1),colormaprange,distplot);
end
