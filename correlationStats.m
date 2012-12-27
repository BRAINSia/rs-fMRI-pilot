[ret, fileStr] = system('ls /paulsen/Experiments/rsFMRI-test/Results/*/*/freesurfer_corr.mat');
index = 1;
fileArray = char();
while (index < length(fileStr))
  [foundStr, count, errmsg, index] = sscanf(fileStr, '%s.mat ');
  fileArray = strvcat(fileArray, foundStr);
  fileStr = fileStr(index+1:end);
end
csvArray = zeros(size(fileArray,1),4);
for index = 1:size(fileArray,1)
    filename = sscanf(fileArray(index, :), '%s');
    load(filename);
    upper = triu(data);
    array = upper;
    for i = 1:numel(array)
        if (array(i) == 0)
            array(i) = nan;
        end
    end
    assert(nanmax(nanmax(array)) == 1);
    % Remove diagonal from statistics
    [rows, columns] = size(array);
    for jj = 1:rows
        array(jj,jj) = NaN;
    end
    minimum = nanmin(nanmin(array));
    maximum = nanmax(nanmax(array));
    median = nanmedian(nanmedian(array));
    mean = nanmean(nanmean(array));
    echo
    newfilename = strcat(filename(1:end-4), '_clean+stats.mat');
    echo
    save(newfilename, 'array', 'minimum', 'maximum', 'median', 'mean');
    % hist(array(:))
    % hist(array(:), 25)
    csvArray(index, 1:4) = [minimum, maximum, median, mean];
    % Write out histogram
    hist(array(:),25)
    title(newfilename)
    print -dpsc2 -append 'corrHist.ps'
    close
end
csvwrite('corrStats.csv', csvArray);