cd 'C:\Users\BNU\Desktop\HBN\Freesurfer'
we=dir(fullfile('hbn*_harmonized.tsv'));
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



