function varargout = iMAP(varargin)
% IMAP MATLAB code for iMAP.fig
%      IMAP, by itself, creates a new IMAP or raises the existing
%      singleton*.
%
%      H = IMAP returns the handle to a new IMAP or the handle to
%      the existing singleton*.
%
%      IMAP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAP.M with the given input arguments.
%
%      IMAP('Property','Value',...) creates a new IMAP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before iMAP_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to iMAP_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help iMAP

% Last Modified by GUIDE v2.5 23-Feb-2015 09:30:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @iMAP_OpeningFcn, ...
    'gui_OutputFcn',  @iMAP_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before iMAP is made visible.
function iMAP_OpeningFcn(hObject, event, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to iMAP (see VARARGIN)

%output logo
if exist('logo_imap.png','file')==2
    axes(handles.axes2);
    try
        imshow('logo_imap.png')
    catch
        imagesc('logo_imap.png')
    end
    
else
    set(handles.axes2,'visible','off');
end
% Choose default command line output for iMAP
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes iMAP wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = iMAP_OutputFcn(hObject, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Preprossesing_data.
function Preprossesing_data_Callback(~, ~, ~)
% hObject    handle to Preprossesing_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%call the preprocessing gui
imap_gui


% --- Executes on button press in Linear_mixed_model.
function Linear_mixed_model_Callback(hObject, eventdata, handles)
% hObject    handle to Linear_mixed_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Asking the user if he wants to see the descriptive statistics
answer_descriptive = questdlg('Would you like to see the descriptive statistics?','','Yes','Skip','Skip');

%% descriptive results
if isempty(answer_descriptive)==0
if strcmp(answer_descriptive,'Yes')==1;

descriptive_part();
else
modelling_part()
end
%% Modelling Analysis
else
clear all
modelling_part()
end

% %%
% clear all
% idx = modelling_analysis();
% if idx ==1 % user did select a file
%     %load  LMMmap
%     [filename, pathname] = uigetfile('*.mat',' Select the FixMap & LMMmap from before to proceed','MultiSelect','on');
%     if iscell(filename)==0 %user selects only multiple file
%         errordlg('FixMap and LMMmap files need to be selected');
%         uiwait(gcf)
%     else
%         fullpathname = strcat(pathname,filename); %user selects only one file
%         if isempty(fullpathname)
%             errordlg('No file selected')
%         else
%             for i = 1:length(filename)
%                 load(strcat(pathname,filename{i}))
%             end
%             if exist('LMMmap','var')&& exist('FixMap','var')
%                 model_plot(LMMmap,FixMap,pathname);
%             else
%                 errordlg('The LMMmap and the FixMap are mandatory to proceed')
%             end
%         end
%     end
%end

% --- Executes on button press in loading_results.
function loading_results_Callback(~, ~, ~)
% hObject    handle to loading_results (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clear all
%answer = questdlg('Would you like to visualize the StatMap or perform modeling,anova, etc');

answer=questdlg('Would you like to visualize the StatMap or perform Linear Contrast?','','StatMap','Linear Contrast','Linear Contrast');
if isempty(answer)==0

if strcmp(answer,'Linear Contrast')==1;
    
    [filename, pathname] = uigetfile('*.mat',' Select FixMap, LMMmap to proceed','MultiSelect','on');
    if iscell(filename)==0 %user selects only multiple file
        h = errordlg('FixMap and LMMmap files need to be selected');
        set(h,'color','white');
        uiwait(gcf)
    else
        fullpathname = strcat(pathname,filename); %user selects only one file
        if isempty(fullpathname)
            h = errordlg('No file selected');
            set(h,'color','white');
        else
            for i = 1:length(filename)
                load(strcat(pathname,filename{i}))
            end
            if exist('LMMmap','var')&& exist('FixMap','var')
                model_plot(LMMmap,FixMap,pathname);
            else
                h = errordlg('The LMMmap and the FixMap are mandatory to proceed');
                set(h,'color','white');
            end
        end
        
    end
else
    [filename, pathname] = uigetfile('*.mat',' Select StatMap to proceed');
    fullpathname = strcat(pathname,filename); %user selects only one file
    
    if isempty(fullpathname)
        errordlg('No file selected')
    else
        load(strcat(pathname,filename))
        if exist('StatMap','var')
            StatMaptoplot=StatMap;
        elseif exist('StatMap_c','var')
            StatMaptoplot=StatMap_c;
        else
            errordlg('The StatMap is mandatory to proceed')
            return
        end
        if strcmp(StatMaptoplot.opt.type,'model')
            %option for LMMdisplay
            %normalized,backgroundfile,colourmap,colormaprange,distplot
            prompt={'background file','colormap','map range for optimal display (0 to 1)','map Value distribution (1 to display)'};
            name=' Options for Display';
            numlines = [1 50;1 50;1 50;1 50];
            defaultanswer={'','''jet''','',''};
            options=inputdlg(prompt,name,numlines,defaultanswer);
            normalized      = 0;
            backgroundfile  = char(options{1});
            cmap            = options{2};
            colormaprange   = str2double(options{3});
            distplot        = str2double(options{4});
        else
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
        end
        %output figure;
        % imapLMMdisplay(StatMap,0)
        if isempty(StatMaptoplot.label)
            StatMaptoplot.label={'unname Contrast'};
        end
        imapLMMdisplay(StatMaptoplot,normalized,backgroundfile(2:end-1),cmap(2:end-1),colormaprange,distplot,[]);
    end
end   
end

% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, ~, ~)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2
%output logo
if exist('logo_imap.png','file')==2
    axes(hObject)
    imshow('logo_imap.png');
else
    set(hObject,'visible','off');
end


% --------------------------------------------------------------------
function uitoggletool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ispc
    winopen('manual_GUI.pdf');
elseif ismac
    system(['open manual_GUI.pdf']);
else
    errordlg('Please open the Guidebook manually in ./matlab/Apps/iMAP/GUI')
end

