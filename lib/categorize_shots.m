function cat = categorize_shots(data,opts)
    header({0,'Categorizing by shot type'})
    cat_labels = unique(data.sync.msr.class);
    num_cats = numel(cat_labels);
    if num_cats ~=0
        cat = cell(num_cats,1);
        for cat_idx = 1:num_cats
            label = cat_labels{cat_idx};
            cat_mask = strcmp(label,data.sync.msr.class);
            cat{cat_idx}.data = struct_mask(data.sync.msr,cat_mask);
        end
%         data.cat.num_cats = num_cats;
    else
        header({2,'No cat data detected.'})
        cat{1}.data = data.sync.msr;
%         data.cat.num_cats = 1;
    end
    header({1,'Done.'})
end