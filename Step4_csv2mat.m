cd '/Users/lizheng/Desktop/同步文件夹/博士研究课题/OHBM会议数据分析/freesurfer_split'
we=dir(fullfile('*_test.csv'));
leng=size(we,1);
for i=1:leng
    loadname = we(i).name;
    data=readtable(loadname, ...
    'FileType', 'text', ...       % 关键：强制按文本文件读取
    'Delimiter', ',', ...        % 明确分隔符为制表符
    'Encoding', 'UTF-8', ...      % 避免中文乱码（可选，需文件编码匹配）
    'ReadVariableNames', true);
    fin = data(:,2:361);
    save([loadname([1:end-4]),'.mat'],'fin')
end



