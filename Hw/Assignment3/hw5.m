a=fliplr(reshape(1:56,7,8));
af_map=mapping(a);
de_map=demapping(af_map);
function map=mapping(a)
map=flip(fliplr(a));
end
function demap=demapping(map)
demap=flip(fliplr(map));
end