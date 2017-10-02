% find waves from ring pressure data 

dt = 1/8;

% 2nd derivative
d2 = [1, -2, 1] / dt.^2;
% third derivative 
d3 = [-1, 2, 0, -2, 1] / dt.^3;




