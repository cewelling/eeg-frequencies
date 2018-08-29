function adjustHist(matrix)
while(1)
    cAxis = xlim;
    xmin = str2double(input('Input desired xmin or press enter if satisfied\n','s'));
    xmax = str2double(input('Input desired xmax or press enter if satisfied\n','s'));
    if isnan(xmin) && isnan(xmax)
        break
%     elseif isnan(xmin)
%         xlim([cAxis(1) xmax])
%     elseif isnan(xmax)
%         xlim([xmin cAxis(2)])
%     else
%         xlim([xmin xmax])
    elseif isnan(xmin)
        hist(matrix(matrix < xmax),100)
    elseif isnan(xmax)
        hist(matrix(matrix > xmin),100)
    end
end
