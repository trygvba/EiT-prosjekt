function [p tri tetr] = loadGeo(filenm)

addpath(genpath('../Geometri'))

p = load([filenm '_nodes.m' ]);
p = p(:,2:4);

tri = load([filenm '_tri.m']);
tri = tri(:,1:3);

tetr = load([filenm '_tetr.m']);
tetr = tetr(:,1:5);

[p tetr tri] = RemoveUnused(p, tetr, tri);

end
