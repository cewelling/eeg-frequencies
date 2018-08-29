for i = 0:9
    files = dir(['stratus105_dartRival*' num2str(i) '.mat']);
    for j = 1:length(files)
        filename = files(j).name
        movefile(filename, [filename(1:37) '_principles.mat'])
    end
end