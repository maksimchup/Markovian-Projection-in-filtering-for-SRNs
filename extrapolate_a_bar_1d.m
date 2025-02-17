function a_bar = extrapolate_a_bar_1d(a_bar, T, varargin)
%fills 'NaN' values in the 'a_bar' based on extrapolation 
% extrapolation on state space, for each fixed time point separately
%
% method = 'linear', 'nearest', ... (see interp1.m)

if ~isempty(varargin)
    method = varargin{1};
else
    method = 'nearest';
end


[r, nx, nt] = size(a_bar);

x_vals = 0:nx-1;
x_vals = reshape(x_vals, [], 1);


for i = 1:r % for each reaction 
    for it = 1:nt
    
        x = x_vals;
    
        v = squeeze(a_bar(i,:,it));
        
        ind = isnan(v);
        x(ind) = [];
        v(ind) = [];
        if length(x) > 1
            a_bar(i,:,it) = interp1(x, v, x_vals, method, 'extrap');
        else
            a_bar(i,:,it) = v;
        end
    end
end

a_bar(isnan(a_bar)) = 0;
a_bar(a_bar < 0) = 0;

end

