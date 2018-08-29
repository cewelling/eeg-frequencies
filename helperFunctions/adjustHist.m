function adjustHist(matrix)
% adjustHist(matrix) takes a matrix of values being used to generate a 
% histogram, prompts the user for x-axis limits, and plots a new histogram
% accordingly. 
% 
% Called from: ftConfigParams.m
% Dependencies: None

while(1)
    % take user input on desired x-axis limits
    xmin = str2double(input('Input desired xmin or press enter if satisfied\n','s'));
    xmax = str2double(input('Input desired xmax or press enter if satisfied\n','s'));
    
    % generate new histogram with data that falls within limits
    if isnan(xmin) && isnan(xmax)
        break
    elseif isnan(xmin)
        matrix = matrix(matrix < xmax);
    elseif isnan(xmax)
        matrix = matrix(matrix > xmin);
    else
        matrix = matrix(matrix > xmin & matrix < xmax);
    end
    hist(matrix,100)
end
