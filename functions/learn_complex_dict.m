function [D, G, err, yy]=learn_complex_dict(x,fn,iter, atoms, sparse) 
%Used to train dictionaries using the K-SVD algorithm. Input signal x and
%dictionary parameters are inputs, and a dictionary, sparse matrix, and
%recovered signals (D*G) are returned.  
%
%Input -  x      : Input signals with individual training signals in the
%                  columns (DEMO data: 10000 X 3001 in the time domain.
%                  3001 training signals of length 10,000 each)  
%         fn     : The range of frequencies that should be used for
%                  training (DEMO data: often 1:1500 for all the
%                  frequencies contained in the 3001 samples)
%         iter   : Number of iterations should be done in training the
%                  dictionary  
%         atoms  : The number of atoms the dictionary should have (Number
%                  of columns in D) 
%         sparse : The number of sparse entries used for the sparse
%                  recovery step (Number of non-zero entries in the columns
%                  of the sparse matrix G. An indicator of how much
%                  sparsity to enforce during learning)
%Output - D      : A dictionary of "size(x, 1)" rows and "num_atom" columns
%                  (DEMO data: dictionary would be 10,000 by num_atom) 
%         G      : The final sparse recovery matrix 
%         err    : All of the RMSE vaules at each iteration 
%         yy     : The input data after frequencies have been removed 

%TAKE FOURIER TRANSFORM
Y = fft(x.'); %Transpose first because fft works on the columns, now have 3001 X 10000 in Fourier domain.
Y(fn(end)+1:end,:) = 0; Y(1:fn(1)-1,:) = 0; yy = real(ifft(Y)); %Remove unwawnted frequencies
Y = Y(1:ceil(end/2),:);

%ASSIGN DICTIONARY LEARNING PARAMETERS
params.data=Y.'; %Transpose back because the input to cksvd should be 10000 X 1501 (after redundant frequencies are removed)
params.memusage='high'; %ksvd parameter, makes computeration faster
params.codemode='sparsity'; %ksvd parameter
params.iternum=iter; %Number of update iterations 
params.Tdata=sparse; %Number of non-zero components for each frequency in the dispersion curve

rndDict=complex(randn(size(params.data,1),atoms),randn(size(params.data,1),atoms)); %Dictates the dictionary size (number of atoms), randomly initialized with complex numbers
params.initdict=rndDict;

%LEARN DICTIONARY
tic %track how long it takes to train
fprintf('\nLet''s train our dictionary.\n'); %tell user what's going on
%Be sure to setup the KSVD algorihtm correctly and includ it in the current directory
[D,G,err] = cksvd(params); %The KSVD algoithm wil train a dictionary by taking the columns of params.data
toc 