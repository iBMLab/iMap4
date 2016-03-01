function [data, subject_in] = function_extension_multiplefile(pathname,filename,extension,handles)
idx = 1; subject_in = 0;
switch extension
    
   case '.txt'

    
 %get information regarding the input file: number of columns and delimitations
[ncolumn,delimitation] = select_delimitation();
h = waitbar(0,'Please wait importing files is in progress','Name','Importing Files','color','w');

switch delimitation
    case 1
 x=[];    
      
 for i = 1:ncolumn
x = strcat(x,'%s');    
 end

readdata = [];

for i = 1:length(filename)
 if ~ishandle(h) %check if the user closes the waitbar
        break
 end
   
readdata1 = textscan(fopen(strcat(pathname,filename{i})),x);
name_file = cell(length(readdata1{1}),1);
name_file(1:length(name_file))= cellstr(filename{i});  
readdata = [readdata;[readdata1,{name_file}]];
waitbar(i/length(filename));
end
if ishandle(h)
delete(h);
end

case 2
 x=[];    
       
for i = 1:ncolumn
x = strcat(x,'%s');    
end
readdata = [];
for i = 1:length(filename)
if ~ishandle(h) %check if the user closes the waitbar
       break
end    
readdata1 = textscan(fopen(strcat(pathname,filename{i})),x,'delimiter','t');
name_file = cell(length(readdata1{1}),1);
name_file(1:length(name_file))= cellstr(filename{i});  
readdata = [readdata;[readdata1,{name_file}]];
waitbar(i/length(filename));
end
if ishandle(h)
delete(h);
end
    
case 3
x=[];    
       
for i = 1:ncolumn
x = strcat(x,'%s');    
end
readdata = [];
for i = 1:length(filename)
if ~ishandle(h) %check if the user closes the waitbar
       break
end    
readdata1 = textscan(fopen(strcat(pathname,filename{i})),x,'delimiter',',');
name_file = cell(length(readdata1{1}),1);
name_file(1:length(name_file))= cellstr(filename{i});  
readdata = [readdata;[readdata1,{name_file}]];
waitbar(i/length(filename));
end
if ishandle(h)
delete(h);
end
    
case 4
x=[];    

for i = 1:ncolumn
x = strcat(x,'%s');    
end
readdata = [];
for i = 1:length(filename)
if ~ishandle(h) %check if the user closes the waitbar
       break
end

readdata1 = textscan(fopen(strcat(pathname,filename{i})),x,'delimiter',';');
name_file = cell(length(readdata1{1}),1);
name_file(1:length(name_file))= cellstr(filename{i});  
readdata = [readdata;[readdata1,{name_file}]];


waitbar(i/length(filename));
end
if ishandle(h)
delete(h);
end
    
    otherwise
        
    msgbox('This delimitation is not treated')
end

for i=1:length(readdata)
e = readdata(:,i);   
AllC = cat(1,e{:});
try
data(:,i) = AllC;    
catch
if ishandle(h)
delete(h);
end
    
errordlg('Please verify the number of columns,the files contain less columns compared to the entered number');
idx = 0;
uiwait(gcf)
end
end

case {'.mat'}
data = []; 
addpath(pathname)
h = waitbar(0,'Please wait importing files is in progress','Name','Importing Files','color','w');

for i = 1:length(filename) 
if ~ishandle(h) %check if the user closes the waitbar
        break
end    
load(strcat(pathname,filename{i}));
variable = fieldnames(load(filename{i}));
if length(variable)==2
errordlg('There is more than one variable saved!');
break;
else
name_file = cell(size(eval(variable{1}),1),1);
name_file(1:length(name_file))= cellstr(filename{i});  
data = [data;[arrayfun(@num2str, eval(variable{1}), 'unif', 0),name_file] ];
waitbar(i /length(filename))
%data = []; % transform data from double to string
% for j = 1:size(data_,2)
% tic
% %data= [data,cellstr(num2str(data_(:,j)))];
% data = [data,data_];
% toc
% end
end
end

if ishandle(h) % if the user did not close the wait bar menu
    delete(h)
end
    
% case {'.xls','.xlsx'}
% for i =1:length(filename)    
% [num,txt,data] = xlsread(strcat(pathname,filename{i}));
% end
% case {'.dat'}
% filetxt=fopen(filename);
% readdata=textscan(filetxt,'%s%s%s%s%s%s%s%s%s%s%s%s%s');
    otherwise
        
data = [];
 errordlg('Sorry this extension is not treated in this version, only matlab (.mat) and text (.txt) extensions are treated.','Data input');        

end


     if idx ==1
    nfiles =size(filename,2);
    answer_user=questdlg(strcat('You selected', {' '}, num2str(nfiles),{' '},'file, Does the Subject column exist?'),'','yes','no','no');
    if strcmp(answer_user,'no')
        
        [data, subject_in] = create_subject(filename,data,handles);
    else
        subject_in = 0;
    end
    
     end    
