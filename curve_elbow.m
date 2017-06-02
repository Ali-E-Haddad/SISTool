function [ind, erat] = curve_elbow(x)
%  [ind, erat] = curve_elbow(x)
%
% Curve elbowing.
% x = length N vector.
% ind  = index of the sample before the elbow.
% erat = ratio of energy in the curve down to ind.

N = length(x);
if N ~= size(x,2)
    x = x';
end
dx       = sum([x(1:N-2);-2*x(2:N-1);x(3:N)]);
[~, ind] = max(dx);
if nargout > 1
    erat = x(1:ind)*x(1:ind)' / (x*x');
end
return