function save_LMM(LMMmap,folder)

prompt = {'Please rename the LMM result',};
dlg_title = 'Saving Variables';
num_lines = 1;
def = {'LMMmap_'};
answer = [];
while isempty(answer)
    answer = inputdlg(prompt,dlg_title,num_lines,def);
end

save(strcat(folder,answer{1}), 'LMMmap','-v7.3');
end
