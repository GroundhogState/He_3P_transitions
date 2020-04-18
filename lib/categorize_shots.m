function cat = categorize_shots(data,opts)
    cli_header(0,'Categorizing by shot type');
    cat_labels = unique(data.check.class);
    num_cats = numel(cat_labels);
    if num_cats ~=0
        cat = cell(num_cats,1);
        for cat_idx = 1:num_cats
            label = cat_labels{cat_idx};
            cat_mask = strcmp(label,data.check.class);
            cat{cat_idx}.data = struct_mask(data.check,cat_mask);
        end
%         data.cat.num_cats = num_cats;
    else
        cli_header(2,'No cat data detected.');
        cat{1}.data = data.check;
%         data.cat.num_cats = 1;
    end
    cli_header(1,'Done.');
end