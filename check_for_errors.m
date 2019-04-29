function data_check = check_for_errors(data,opts)
    

    low_count_mask = data.cal.mdl > opts.check.min_counts;
%     data.cal.probe_set = data.sync.msr.probe_set(low_count_mask);
%     data.cal.signal = data.cal.signal_premask(low_count_mask);


    on_mask = data.sync.msr.exposure> opts.ai.exposure;
    power_mask = data.sync.msr.mean_power > opts.ai.pd_setpoint;
    ai_mask = on_mask & power_mask;
    if opts.check.ai_override
        ai_mask = ones(size(ai_mask));
    end
    fprintf('PD performance failure in %u of %u shots\n',sum(~ai_mask),length(ai_mask));

    master_mask = data.cal.low_count_mask & ai_mask;

    data_check = struct_mask(data.cal,master_mask);

    fprintf('Ignoring %u shots with atom number<%u\n',sum(~data.cal.low_count_mask),opts.check.min_counts)
    fprintf('Ignoring %u shots with probe beam failure\n',sum(~ai_mask))
    fprintf('Ignoring %u shots in total\n',sum(~master_mask))
end