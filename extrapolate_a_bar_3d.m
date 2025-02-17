function a_bar = extrapolate_a_bar_3d(a_bar, T, varargin)
%fills 'NaN' values in the 'a_bar' 
% extrapolation on space for each time point

if ~isempty(varargin)
    method = varargin{1};
else
    method = 'nearest';
end

[r, nx, ny, nz, nt] = size(a_bar);

x_vals = 0:nx-1;
y_vals = 0:ny-1;
z_vals = 0:nz-1;
[z_vals,y_vals,x_vals] = meshgrid(z_vals,y_vals,x_vals);

z_vals = reshape(z_vals, [], 1);
x_vals = reshape(x_vals, [], 1);
y_vals = reshape(y_vals, [], 1);
    

for i = 1:r % for each reaction 
    for it = 1:nt % for each time point
        x = x_vals;
        y = y_vals;
        z = z_vals;
    
        v = squeeze(a_bar(i,:,:,:,it));
        v = reshape(v, [], 1);
        
        ind = isnan(v);
        x(ind) = [];
        y(ind) = [];
        z(ind) = [];
        v(ind) = [];
        try 
            F = griddedInterpolant(x,y,z,v,method,method);
            res = F(x_vals,y_vals,z_vals);
            a_bar(i,:,:,:,it) = reshape(res, nx, ny, nz);
        catch
        end
    end
    a_bar(i,isnan(a_bar(i,:))) = nanmean(a_bar(i, :));
end


a_bar(a_bar < 0) = 0;

end

