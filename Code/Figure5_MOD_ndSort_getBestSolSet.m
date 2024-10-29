clear all
clc

% 更改工作路径
new_filePath = '******';
cd(new_filePath);

% 显示新的工作路径
new_current_filePath = pwd;
disp(['新的工作路径是：', new_current_filePath]);

% 首先，指定包含csv文件的文件夹
filePath = strcat(new_filePath,'\opt30ResMOD');
% 然后，使用dir函数列出文件夹中所有csv文件的名称
fileNameList = dir(fullfile(filePath, '*.csv'));

for i = 1:length(fileNameList)
    % 获取当前文件的名称
    currentFileName = fileNameList(i).name;
    rsm = strsplit(currentFileName,'_');
    rsm = string(rsm(2));
    
    % 使用readtable函数将当前文件导入到MATLAB表格中
    currentPath = fullfile(filePath, currentFileName);
    all_pf = readtable(currentPath);
    
    all_pf(1,:) = []; % 删除第一行和第一列
    all_pf(:,1) = [];
    all_pf = table2array(all_pf);   
    
    O = all_pf(:,318:319);
    V = all_pf(:,1:317);
    V = floor(V);
    
    all_pf = [V O];
    unique_all_pf=unique(all_pf, 'rows');
    f=non_domination_sort_mod(unique_all_pf,2,317);
    best_pf=f(f(:,320)==1,:);
    
    subplot(2,4,i);
    scatter(unique_all_pf(:,318),-unique_all_pf(:,319),15,'b');
    hold on;
    scatter(best_pf(:,318),-best_pf(:,319),15,'r','filled');

    ylabel(rsm+'(-)');
    xlabel('Cost(M€)');
        
    writeFileName = strcat(new_filePath,'\ndSortResMOD\MOD_',rsm,'_ndSortRes.xlsx');
	sheetName = 'best_pf';
    % writematrix(all_pf, writeFileName, 'Sheet', 'all_pf');
    % writematrix(unique_all_pf, writeFileName, 'Sheet', 'unique_all_pf');
    % writematrix(f, writeFileName, 'Sheet', 'f');
    writematrix(best_pf, writeFileName, 'Sheet', sheetName);

end
