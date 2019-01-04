function V = omp( A, X, K ) 
%OMP  Computes the sparse orthogonal matching pursuit solution
%   V = OMP(A, B, K)  computes the sparse solution V given a basis matrix 
%   A and observations X.
%
%   INPUTS: 
%       A: An M-by-N basis matrix
%       X: An M-by-Q vector of observations
%       K: Number of sparse components
%
%  OUTPUTS: 
%       V: The N-by-Q wavenumber basis pursuit denoising results
%

    % DEFINE SIZES 
    N   = size(A,2);            % Number of atoms
    M   = size(A,1);            % Aize of atoms
    Q   = size(X,2);            % Number of outputs/frequencies

    % INITIALIZE VARIABLES
    V         = zeros(N,Q);     % Solution
    V_T       = zeros(K,Q);     % Solution with only non-zero indices
    indx_set  = zeros(K,Q);     % Indices with non-zero values
    atoms     = zeros(M,K,Q);   % Chosen dictionary atoms for each frequency

    % INIATIVE ALGORITHM 
    r   = X;                    % Initial residual
    Vr  = A'*r;                 % Initial solution from residual xr
    
    % LOOP OVER NUMBER OF SPARSE COMPONENTS
    for k = 1:K
        
        % FIND CORRESPONDING INDICES
        [~,ind_new]   = max(abs(Vr), [], 1);             % Find best match
        indx_set(k,:) = ind_new;                         % Add to index set
        atoms(:,k,:)  = permute(A(:,ind_new), [1 3 2]);  % Get cooresponding atom

        % UPDATE RESIDUAL
        for q = 1:Q  % Loop over outputs
            V_T(1:k,q) = atoms(:,1:k,q) \ X(:,q);          % Find least-squares fit
            V( indx_set(1:k,q), q )   = V_T(1:k,q);        % Places results in full vector
            r(:,q)  = X(:,q) - atoms(:,1:k,q)*V_T(1:k,q);  % Find new residual
        end
        
        % COMPUTE SOLUTION FROM RESIDUAL Vr FOR NEXT ITERATION
        if k < K, Vr  = A'*r; end

    end

end
