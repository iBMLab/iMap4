function [ncolumn,delimitation] = select_delimitation()

prompt = {'Enter the number of column:      ';'Select the delimitation of the file:'};
name = '';
formats = struct('type', {}, 'style', {}, 'items', {},'format', {}, 'limits', {}, 'size', {});
formats(1,1).type   = 'edit';
formats(1,1).format = 'integer';
formats(1,1).size = 70;
formats(2,1).type   = 'list';
formats(2,1).style  = 'popupmenu';
formats(2,1).items  = {'Space','Tabulation','Comma','Semicolon'};
defaultanswer = {12, 1};
formats(2,1).size = 70;
answer = inputsdlg(prompt, name, formats, defaultanswer);
ncolumn = answer{1}; delimitation = answer{2};
