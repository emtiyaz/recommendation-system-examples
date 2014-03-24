function y = apply_quadratic(Q,x)

[n,d] = size(x);
y = zeros(n,1);
for i=1:n
    y(i) = x(i,:)*Q{1}*x(i,:)' + x(i,:)*Q{2}(:) + Q{3};
end


% C:\Documents and Settings\gbouchar\My Documents\work\mytools\statlearn\utils\mathfun