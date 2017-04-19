function [ind, eperc] = curve_elbow(x)
%  [ind, eperc] = curve_elbow(x)
%
% Curve elbowing.
% x = length N vector.
% ind   = index of the sample before the elbow.
% eperc = ratio of energy in the curve down to ind.
%
% Coded by: Ali Haddad

N = length(x);
if N ~= size(x,2)
    x = x';
end
dx       = sum([x(1:N-2);-2*x(2:N-1);x(3:N)]);
[~, ind] = max(dx);
eperc    = x(1:ind)*x(1:ind)' / (x*x');
return