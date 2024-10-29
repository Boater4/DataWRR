clear all
clc

% 更改工作路径
new_filePath = '******';
cd(new_filePath);

% 显示新的工作路径
new_current_filePath = pwd;
disp(['新的工作路径是：', new_current_filePath]);

% 首先，指定包含csv文件的文件夹
filePath = strcat(new_filePath,'\opt30ResHAN');
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
    
	O = all_pf(:,35:36);
	V = floor(all_pf(:,1:34));
    
    all_pf = [V O];
    unique_all_pf=unique(all_pf, 'rows');
	f=non_domination_sort_mod(unique_all_pf,2,34);
	best_pf=f(f(:,37)==1,:);
	
    subplot(6,6,i);	
	scatter(unique_all_pf(:,35),-unique_all_pf(:,36),15,'b');
	hold on;
	scatter(best_pf(:,35),-best_pf(:,36),15,'r','filled');
    ylabel(rsm+'(-)');
    xlabel('Cost(M$)');
    
    str1 = '\ndSortResHAN\HAN';
    str2 = '_';
    str3 = '_ndSortRes.xlsx';
    writeFileName = strcat(new_filePath,str1,str2,rsm,str3); %,string(i-1)
	sheetName = 'best_pf';
    % writematrix(all_pf, writeFileName, 'Sheet', 'all_pf');
    % writematrix(unique_all_pf, writeFileName, 'Sheet', 'unique_all_pf');
    % writematrix(f, writeFileName, 'Sheet', 'f');
    writematrix(best_pf, writeFileName, 'Sheet', sheetName);

end
