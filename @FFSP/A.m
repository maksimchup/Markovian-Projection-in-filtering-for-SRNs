function A = A(filter, t, Yt)
%Creates matrix 'A' for the filtering equation at time point 't'

% indices & values of elements in sparse matrix 'A'
i_arr = ones(1, filter.N*filter.model.r);
j_arr = ones(1, filter.N*filter.model.r);
v_arr = zeros(1, filter.N*filter.model.r);

k = 1; % ind of current element

Xi = zeros(filter.model.d-numel(filter.observed_ind), 1);

for i = 1:filter.N % row index

    % full state vector
    Zi = Xi;
    for io = 1:numel(Yt)
        Zi = insert(Zi, filter.observed_ind(io), Yt(io));
    end
    ai = filter.model.a(Zi, t);
    
    % diagonal element
    i_arr(k) = i;
    j_arr(k) = i;
    v_arr(k) = -sum(ai);
    k = k + 1;


    % off-diagonal elements
    for ir = filter.unobs_reac_ind 
        Zj = Zi - filter.model.V(:, ir);

        Xj = remove(Zj, filter.observed_ind);       
        j = filter.state2ind(Xj);

        if j ~= 0 && j ~= i
            aj = filter.model.a(Zj, t);
            
            i_arr(k) = i;
            j_arr(k) = j;
            v_arr(k) = aj(ir);
            k = k + 1;
        end
    end

    Xi = filter.next_state(Xi);
end

A = sparse(i_arr, j_arr, v_arr);

end


function new_arr = insert(arr, i, value)
%insert value in array (size (~,1))
new_arr = [arr(1:i-1); value; arr(i:end)];
end

function new_arr = remove(arr, i)
%remove element from array (size (~,1))
new_arr = arr;
new_arr(i) = [];
end