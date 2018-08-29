% peak_find_test_main
q_factor = 5; % for notch filter
pk_threshold_set = [0.2E-4]; % local pk amplitude
pk_threshold_local_set = [0.3E-4]; % local pk amplitude

Filter_switch = 1;
FBS = [0.8];
y0 = HRnotch_filter(y0,q_factor);
if Filter_switch==1
    %     y0 = BSfilter(y0',4,FBS,[],'but','twopass','no');
    y0 = LPfilter(y0',4,FBS,[],'but','twopass','no');
end
[locs, peakMag] = peakfinder(y0,pk_threshold_local_set,-pk_threshold_set,-1,0);
