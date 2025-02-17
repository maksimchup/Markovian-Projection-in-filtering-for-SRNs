function ind = state2ind(filter, X)
%index of state 'X' in vector 'pi'

if any(X > filter.max_X) || any(X < 0)
    ind = 0;
    return 
end

ind = filter.ind_offsets*X + 1;

end