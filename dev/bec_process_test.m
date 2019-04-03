% bec_transition_test_fun

%Create a synthetic data set with well known peak structure and verify that
%the function returns the correct values!
opts = [];
data.lv = [];
data.wm = [];
data.tdc = [];

test_fun = bec_transition_process(opts,data);
test_fun == test_values?