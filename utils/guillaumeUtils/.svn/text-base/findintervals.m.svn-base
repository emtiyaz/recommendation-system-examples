function pushes = findintervals(points,bounds)
% findintervals - find the index of the intervals
% I = findintervals(points,bounds)
% 
% %example:
% findintervals([.1 .5 1.1 1.2 2 2.5 3.5],[1 2 3])
%
% example 2:
% a = [.1 .5 1.1 1.2 1.7 2.5 3.5]; b = [1 2 3 4];
% i = findintervals(a,b);min(b(i)-a);
% plot(a,1:length(a));line(b(i),1:length(a),'Color','r');
% 
% WARNING: bounds must be sorted
% WARNING: points must not be nan



np = length(bounds);
ne = length(points);
[times,lbl] = cell2lbl({points(:),bounds(:)+eps});

[u,ord] = sort(times);

loc = find(ord>ne)-(1:np)';%date of bounds
a = full(sparse(ones(np,1),loc+1,ones(np,1),1,ne+1));
pushes = cumsum(a)+1;
pushes = pushes(1:end-1);



