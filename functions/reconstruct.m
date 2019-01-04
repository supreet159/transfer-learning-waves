function [ V1,xx ] = reconstruct(test_data, full_dictionary, indxs, fn, spars)
%Reconstructs a full wavefield "xx" from a partial sampling contained in
%"test_data", using orthagonal matching pursuit (OMP).
%
%Input  - test_data: Data from experiment or data after sampling [time X space]
%         full_dictionary: Dictionary used for reconstruction {space X no. of atoms}
%         indxs: Indices of the the locations that the samples in test_data are drawn from.
%         fn: Frequancies used in the reconstruction
%         spars: Sparsity for OMP
%         num_atoms: Number of atoms in the dictionary "D"
%         Q: Number of time samples
%Output - xx: Reconstructed data (DEMO data: 10,000 by 3001)
%         V1: Sparse coefficiencts, sparse matrix


% DEFINE CONSTANTS
Q=size(test_data,1);

% COMPUTE FOURIER TRANSFORM
Xtest = complex(fft(test_data));

% NORMALIZE DICTIONARY
newdictionary = full_dictionary(indxs,:);
nf = 1./sqrt(sum(abs(newdictionary).^2,1)); % normalize the dictionary
newdictionaryN = bsxfun(@times, newdictionary, nf);

% BUILD GRAMMIAN
Gm = newdictionaryN'*newdictionaryN;

% RECOVER SPARSE REPRESENTATION
V1 = comp(newdictionaryN, (Xtest(fn,:)).', Gm, spars);

% RECONSTRUCT DATA 
Z = zeros(Q,size(newdictionary,1));
Z(fn,:) = (newdictionary*bsxfun(@times, V1, nf(:))).';  
xx = 2*real(ifft(Z));

end


function [X] = dictionary_sws(dictionary, SparseMatrix )
   X = (dictionary*SparseMatrix).';   
end