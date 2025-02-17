function x = next_state(filter, x)
%get next state in vector 'pi'
% Let 'x' have index 'i', this funcion returns the state corresponding to 
% index 'i+1'

x(1) = x(1) + 1;
for i = 1:length(x)-1
    if x(i) > filter.max_X(i)
        x(i) = 0;
        x(i+1) = x(i+1) + 1;
    end
end

end