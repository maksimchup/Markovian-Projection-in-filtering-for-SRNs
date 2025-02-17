function filter = resample(filter)
%Resampling step in PF
% ref: https://doi.org/10.1007/978-0-387-76896-0   (section 9.2.1)
% -------------------------------------------------------------------------
%INPUT
% filter      : object of class PF
% -------------------------------------------------------------------------
%OUTPUT
% filter      : PF with updated states & weights
% -------------------------------------------------------------------------



Z = filter.Z;
w = filter.w;
N = filter.M;

w = w / sum(w);

u = rand(N-1, 1);

g = N;
h = N;
o = zeros(N, 1);

for i = 1:N-1
    Nw = N*w(i);
    Nwint = floor(Nw);
    Nwfrac = Nw - Nwint;
    gint= floor(g);
    gfrac = g - gint;

    if Nwfrac + (g-Nw) - floor(g-Nw) < 1
        if u(i) < 1 - Nwfrac/gfrac
            o(i) = Nwint;
        else
            o(i) = Nwint + h - gint;
        end
    else
        if u(i) < 1 - (1-Nwfrac)/(1-gfrac)
            o(i) = Nwint + 1;
        else
            o(i) = Nwint + h - gint;
        end
    end
    g = g - Nw;
    h = h - o(i);
end

o(N) = h;


j = 1;
for i = 1:N
    filter.w(i) = 1;
    for l = 1:o(i)
        filter.Z(:, j) = Z(:, i);
        j = j + 1;
    end
end


end