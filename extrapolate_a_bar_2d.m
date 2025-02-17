function a_bar = extrapolate_a_bar_2d(a_bar, T, varargin)
%fills 'NaN' values in the 'a_bar' 
% extrapolation on space (1d) + time (1d) simultaneously 

if ~isempty(varargin)
    method = varargin{1};
else
    method = 'nearest';
end

[r, nx, ny, nt] = size(a_bar);

x_vals = 0:nx-1;
y_vals = 0:ny-1;
t_vals = linspace(0, T, nt);

[y_vals,x_vals,t_vals] = meshgrid(y_vals,x_vals,t_vals);

x_vals = reshape(x_vals, [], 1);
y_vals = reshape(y_vals, [], 1);
t_vals = reshape(t_vals, [], 1);
    

for i = 1:r % for each reaction 
    x = x_vals;
    y = y_vals;
    t = t_vals;

    v = squeeze(a_bar(i,:,:,:));
    v = reshape(v, [], 1);
    
    ind = isnan(v);
    x(ind) = [];
    y(ind) = [];
    t(ind) = [];
    v(ind) = [];
    
    F = scatteredInterpolant(x,y,t,v, method);
    res = F(x_vals, y_vals, t_vals);
        
    a_bar(i,:,:,:) = reshape(res, nx, ny, nt);
end

a_bar(isnan(a_bar)) = 0;
a_bar(a_bar < 0) = 0;

end

