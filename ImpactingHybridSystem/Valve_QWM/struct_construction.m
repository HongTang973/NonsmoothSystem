function  S = struct_construction(keys, values)

n = max(size(keys,1),size(keys,2));
S = [];
for i  =1 : n
    S.(string(keys(i))) = values(i);
end
end