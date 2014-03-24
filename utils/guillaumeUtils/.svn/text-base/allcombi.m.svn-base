function C = allcombi(sz)
%ALLCOMBI
%
%   Generate all the possible combinations of numbers between 1 and sz(i) for each column i
%   The vector C of row generated is therefore of size prod(sz)
%   Example : 
%   allcombi([2 3]) gives
%   1 1
%   1 2
%   1 3
%   2 1
%   2 2
%   2 3
%
%   Example 2:
%   L=[2 3 4];(1+(allcombi(L)-1)*[1 cumprod(L(1:end-1))]')'

dim = length(sz);
% if dim>100
%     C = ones(sum(sz)-dim+1,length(sz));
%     pos = 2;
%     for i=1:length(sz)        
%         C(pos:pos+sz(i)-2,i)=(2:sz(i))';        
%         pos = pos+sz(i)-1;
%     end
%     return
% end
switch(dim)
case 1
    C = (1:sz)';
case 2
   [r1,r2] = meshgrid(1:sz(1),1:sz(2));
   C = [r1(:) r2(:)];
otherwise
    C = zeros(prod(sz),length(sz));
    for i=1:length(sz)
       u=[1:sz(i)];
       v=(u'*ones(1,prod(sz(i+1:end))))';
       w = v(:)*ones(1,prod(sz(1:i-1)));
       C(:,i) = w(:);
    end
end


