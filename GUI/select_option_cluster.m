function mccopt  = select_option_cluster(mccopt)

quest1 = 'Input the cluster size threshold, only the Cluster size bigger than the threshold will be considered as significant';
quest2 = 'Input the cluster mass threshold, only the Cluster mass (sum of statistic value inside the cluster) larger than the threshold will be considered as significant';
prompt = {quest1,quest2};
dlg_title = 'Cluster';
num_lines = 1;
def = {'',''};
answer = inputdlg(prompt,dlg_title,num_lines,def);

mccopt.clustSize = str2double(answer{1});
mccopt.clustVal  = str2double(answer{2});

if isempty(answer{1})==0 && isempty(answer{2}) ==0

warndlg('iMAP 4 will select the cluster satisfy both criteria')
end