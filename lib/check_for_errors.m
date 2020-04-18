function data_check = check_for_errors(data,opts)
    cli_header(0,'Performing error checks');

    low_count_mask = data.cal.mdl > opts.check.min_counts;
%     data.cal.probe_set = data.sync.msr.probe_set(low_count_mask);
%     data.cal.signal = data.cal.signal_premask(low_count_mask);

    on_mask = data.sync.msr.exposure> opts.ai.exposure;
    power_mask = data.sync.msr.mean_power > opts.ai.pd_setpoint;
    ai_mask = on_mask & power_mask;
    if isfield(data.ai,'error')
        warning("No analog in data recorded! Automatically overriding, don't trust data.")
        opts.check.ai_override = true;
    end
    if isfield(opts.check,'ai_override')
        if opts.check.ai_override
            ai_mask = ones(size(ai_mask));
        end
    end
    fprintf('PD uptime below %.4f in %u of %u shots\n',opts.ai.exposure,sum(~on_mask),length(on_mask));
    fprintf('PD setpt below %.4f %u of %u shots\n',opts.ai.pd_setpoint,sum(~power_mask),length(power_mask));
    fprintf('Nominal diagnostic laser power %.3f +- %.3f mw\n',1e3*nanmean(data.sync.msr.mean_power*opts.ai.volt_to_watt_factor),1e3*nanstd(data.sync.msr.mean_power*opts.ai.volt_to_watt_factor));
    master_mask = data.cal.low_count_mask & ai_mask;

    data_check = struct_mask(data.cal,master_mask);

    fprintf('Ignoring %u shots with atom number < %u\n',sum(~data.cal.low_count_mask),opts.check.min_counts)
%     fprintf('Ignoring %u shots with probe beam failure\n',sum(~ai_mask))
    fprintf('Ignoring %u shots in total\n',sum(~master_mask))
    cli_header(1,'Done.');
end