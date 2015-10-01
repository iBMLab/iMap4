function[ data,subject_in] = function_extension_file(pathname,filename,extension)

switch extension
    
    case '.txt'
        
        filetxt=fopen(strcat(pathname,filename));
        
        %get information regarding the input file: number of columns and delimitations
        [ncolumn,delimitation] = select_delimitation();
        
        switch delimitation
            case 1
                x=[];
                
                for i = 1:ncolumn
                    x = strcat(x,'%s');
                end
                readdata = [];
                readdata1=textscan(filetxt,x);
                name_file = cell(length(readdata1{1}),1);
                name_file(1:length(name_file))= cellstr(filename);
                readdata = [readdata;[readdata1,{name_file}]];
                
            case 2
                x=[];
                
                for i = 1:ncolumn
                    x = strcat(x,'%s');
                end
                readdata = [];
                readdata1=textscan(filetxt,x,'delimiter','t');
                name_file = cell(length(readdata1{1}),1);
                name_file(1:length(name_file))= cellstr(filename);
                readdata = [readdata;[readdata1,{name_file}]];
            case 3
                x=[];
                
                for i = 1:ncolumn
                    x = strcat(x,'%s');
                end
                readdata = [];
                readdata1=textscan(filetxt,x,'delimiter',',');
                name_file = cell(length(readdata1{1}),1);
                name_file(1:length(name_file))= cellstr(filename);
                readdata = [readdata;[readdata1,{name_file}]];
            case 4
                x=[];
                
                for i = 1:ncolumn
                    x = strcat(x,'%s');
                end
                readdata = [];
                readdata1=textscan(filetxt,x,'delimiter',';');
                name_file = cell(length(readdata1{1}),1);
                name_file(1:length(name_file))= cellstr(filename);
                readdata = [readdata;[readdata1,{name_file}]];
        end
        
        % concatenate the data
        data= cat(2,readdata{:});
        addpath(pathname)
    case {'.mat'}
        load(strcat(pathname,filename))
        variable = fieldnames(load(strcat(pathname,filename)));
        if length(variable)==2
            errordlg('There is more than one variable saved!');
        else
            % data = []; % transform data from double to string
            % data_ = eval(variable{1});
            % for j = 1:size(data_,2)
            % data= [data,cellstr(num2str(data_(:,j)))];
            % end
            data = [];
            name_file = cell(length(eval(variable{1})),1);
            name_file(1:length(name_file))= cellstr(filename);
            data = [data;[arrayfun(@num2str, eval(variable{1}), 'unif', 0),name_file] ];
        end
        
        % case {'.xls','.xlsx'}
        % [num,txt,data] = xlsread(strcat(pathname,filename));
        
        % case {'.dat'}
        % filetxt=fopen(filename);
        % readdata=textscan(filetxt,'%s%s%s%s%s%s%s%s%s%s%s%s%s');
        
    otherwise
        data = [];
        errordlg('Sorry this extension is not supported in this version, only matlab (.mat) and text (.txt) files are supported.','Data input');
        
end



nfiles =size(filename,1);
answer_user=questdlg(strcat('You selected', {' '}, num2str(nfiles),{' '},'file, Does the Subject column exist?'),'','yes','no','no');
if strcmp(answer_user,'no')
    
    errordlg('Unable to proceed, the File does not include subject predictor','Subject');
    
else
    subject_in = 0;
end