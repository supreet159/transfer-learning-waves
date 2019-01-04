function [ D, V ] = ksvd_simple( X, D, K, iter )
%KSVD_SIMPLE Summary of this function goes here
%   Detailed explanation goes here

% COMPUTE CONSTANT LENGTH
M = size(X, 1);
N = size(D, 2); 
Q = size(X, 2);

% NORMALIZE DICTIONARY
D = D./sqrt(sum(abs(D).^2, 1));

for kk = 1:iter
    
    % PERFORM OMP WITH A SINGLE SPARSITY
    V = omp( D, X, K );

    % DETERMINE THE ATOMS USED AND UNUSED
    atoms_used = double(abs(sum(V,2))~=0).*(1:N).';
    atoms_used = atoms_used(atoms_used ~= 0);
    atoms_unused = double(abs(sum(V,2))==0).*(1:N).';
    atoms_unused = atoms_unused(atoms_unused ~= 0);

    % KEEP TRACK OF UNUSED TRAINING DATA
    sigs_unused = 1:Q;

    % GET RANDOM PERMUTATION
    p1 = randperm(length(atoms_used));
    p2 = randperm(length(atoms_unused));

    % DICTIONARY LEARNING
    for j = 1:length(p1)

        % DETERMINE INDEX OF ATOM BEING ANLAYZED
        indx = atoms_used(p1(j));

        % EXTRACT DICTIONARY ATOM AND SPARSE REPRESENTATION
        Dj = D(:,indx);
        Vj = V(indx,:);

        % TRUNCATE SPARSE REPRESENTATION AND DATA
        sig_indx = find(sum(abs(Vj),1)~=0);
        Vj = Vj(:,sig_indx);
        Vs = V(:,sig_indx);
        Xs = X(:,sig_indx);

        % COMPUTE ERROR AND ESTIMATE ATOM
        E = Xs - D*Vs + Dj*Vj;
        %[atom, Sat, Vat] = svds(E,1);
        atom = E*Vj';
        atom = atom/norm(atom);

        % ESTIMATE NEW SPARSE VALUE
        %Vju = omp( E, atom, 1 );
        Vju = atom'*E;
        [~, ii] = max(abs(Vju));
        sigs_unused = setdiff(sigs_unused, sig_indx(ii));


        D(:,indx) = atom;
        V(indx,sig_indx) = Vju;

    end

    % CLEAN-UP 
    for j = 1:length(p2)

        indx = atoms_unused(p2(j));

        E = sum(abs(X(:,sigs_unused) - D*V(:,sigs_unused)).^2);
        [d, i] = max(E);
        atom = X(:,sigs_unused(i));
        atom = atom/norm(atom);
        Vju = 0;

        sig_indx = sigs_unused(i);
        sigs_unused = setdiff(sigs_unused, sigs_unused(i));

        D(:,indx) = atom;
        V(indx,sig_indx) = Vju;


    end

    
    V = omp( D, X, K );
    Xhat = D*V;
    cc = real(sum(sum(X.*conj(Xhat))))./norm(X, 'fro')/norm(Xhat, 'fro');
    
    fprintf('ITERATION %i: Correlation Coefficient: %f (%i unused atoms)\n', kk, cc, length(atoms_unused))
    
end

end


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
