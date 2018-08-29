matchColumn = [];

for i = 1:size(A_lowTrans)
    if ~isnan(A_lowTrans(i,1)) && A_lowTrans(i,1) ~= 0
        [a, b] = find(A_lowSegs == A_lowTrans(i,1));
        if isempty(a)
            disp(['::: ' num2str(i) ' :::'])
        else
            matchColumn = [matchColumn b];
        end
    end
end