function out = struct_mask(s,m)
% Assuming s is a structure with equal sizes in each field, mask the
% fields according to m

fieldnames = fields(s);
nfields = numel(fieldnames);
out = [];
for ii=1:nfields
   this_field = fieldnames{ii};
   this_data = getfield(s,this_field);
%    this_data = s(this_field);
   out = setfield(out,this_field,this_data(m));
end


end