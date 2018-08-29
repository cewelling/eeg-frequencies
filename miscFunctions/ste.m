function error = ste(data)

nData = sum(~isnan(data));
error = nanstd(data)./sqrt(nData-1);

if length(data) < size(data,2)
    sprintf('Hold on! Are these data the right dimensions for ste?')
end

