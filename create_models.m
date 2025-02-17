%This script creates .mat files for SRN models 


%%  bistable gene network model (8d)
% species: [G1, G2, G1*, G2*, mRNA1, mRNA2, Protein1, Protein2]
% ref: section III.B at https://doi.org/10.1063/1.5021242

% stoichiometric matrix
V = [   
    [0; 0; 0; 0; +1; 0; 0; 0],   ...  0 -> mRNA1
    [0; 0; 0; 0; 0; +1; 0; 0],   ...  0 -> mRNA2
    [0; 0; 0; 0; -1; 0; 0; 0],   ...  mRNA1 -> 0
    [0; 0; 0; 0; 0; -1; 0; 0],   ...  mRNA2 -> 0
    [0; 0; 0; 0; 0; 0; +1; 0],   ...  mRNA1 -> mRNA1 + Protein1
    [0; 0; 0; 0; 0; 0; 0; +1],   ...  mRNA2 -> mRNA2 + Protein2
    [0; 0; 0; 0; 0; 0; -1; 0],   ...  Protein1 -> 0
    [0; 0; 0; 0; 0; 0; 0; -1],   ...  Protein2 -> 0
    [0; +1; 0; -1; 0; 0; 0; 0],  ...  Protein1 + G2* -> Protein1 + G2
    [+1; 0; -1; 0; 0; 0; 0; 0],  ...  Protein2 + G1* -> Protein2 + G1
    [0; 0; 0; 0; +1; 0; 0; 0],   ...  G1* -> G1* + mRNA1
    [0; 0; 0; 0; 0; +1; 0; 0],   ...  G2* -> G2* + mRNA2
    [-1; 0; +1; 0; 0; 0; 0; 0],  ...  G1 -> G1*
    [0; -1; 0; +1; 0; 0; 0; 0],  ...  G2 -> G2*
    [+1; 0; -1; 0; 0; 0; 0; 0],  ...  G1* -> G1
    [0; +1; 0; -1; 0; 0; 0; 0],  ...  G2* -> G2
    ];

% reaction rates
theta = [0.1; 0.1; 0.05; 0.05; 5; 5; 0.2; 0.2; ...
        0.1; 0.1; 1; 1; 0.03; 0.03; 1e-6; 1e-6];

% propensities
a = @(z,t) [ 
    theta(1);
    theta(2);
    theta(3)*z(5);
    theta(4)*z(6);
    theta(5)*z(5);
    theta(6)*z(6);
    theta(7)*z(7);
    theta(8)*z(8);
    theta(9)*z(7)*z(4);
    theta(10)*z(8)*z(3);
    theta(11)*z(3);
    theta(12)*z(4);
    theta(13)*z(1);
    theta(14)*z(2);
    theta(15)*z(3);
    theta(16)*z(4);
    ];

model = SRN(V, a);
model.is_homogeneous = true;

% Initial state
Z0 = [0; 0; 1; 1; 0; 0; 0; 0];

filename = 'Data/models/bistable_gene_network.mat';
save(filename, 'model', 'Z0');




%% Linear signalling cascade models (d=3..8)
% species: [Z_1, ..., Z_d]
% ref: https://doi.org/10.1371/journal.pcbi.1009623

for d = 3:8 % # of species
    
    % stoichiometric matrix
    V = [
        eye(d) - diag(ones(d-1, 1), 1)   ... Z_i -> Z_{i+1} (including 0 -> Z_1)
        -eye(d),                         ... Z_i -> 0
        ];
    
    % reaction rates
    theta = [10; 5*ones(d-1, 1); 1*ones(d, 1)];
    
    % propensities
    a = @(z,t) [ 
        theta(1);
        theta(2:d).*z(1:d-1);
        theta(d+1:end).*z;
        ];
    
    model = SRN(V, a);
    model.is_homogeneous = true;
    
    % Initial state
    Z0 = [5; 5; 4; zeros(d-3, 1)];
    
    filename = ['Data/models/linear_cascade_' num2str(d) '.mat'];
    save(filename, 'model', 'Z0');

end



