function corrs = recovery_correlation(D, data_set, fn, num_sparse, time_cut, cortest, scltest)
%Creates recover correlations by taking each sampling in "data_set",
%perfomring a recovery using the dictionary "D" and the "reconst" function.
%The corresponding indices from the sampling are then correlated with the
%recovery. One correlation is produced from each sampling. These are them
%recoupled with the label indicating damage state and returned as
%"outcome". 

%Input:  D: The dictionary used for the reconstruction
%        data_set: The set of samplings and labels produced using "create_data_set"
%        fn: The frquencies that will be used in the reconstruction
%        num_sparse: The number of sparse elements used in the reconstruction
%        time_cut: The number of time samples that will actually be used for the correlation
%Output: outcome: List of true labels coupled with correlations between the sampling and corresponding indices of the recovered sampling. 

fprintf('Recovery correlations: '); % Tell user what is going on

%  INITIALIZE VARIBALES
C = cortest+scltest;
corrs = zeros(data_set.num_samplings,2+C); % Label and correlation

% LOOP OVER NUMBER OF SAMPLINGS
fprintf(repmat(' ', 1, 22+8*C));
for i = 1:size(corrs,1)
    
    % COLLECT DATA
    data = data_set.(['s',num2str(i)]);
    train_data = data_set.(['t',num2str(i)]);
    random_rows = data_set.(['i',num2str(i)]);
    corrs(i,1) = data_set.(['l',num2str(i)]);
    
    % RECONSTRUCT DATA
    [~,xx] = reconstruct(data.',D,random_rows,fn,num_sparse); xx=xx.'; %reconstruct sampling
    h = xx(:,1:time_cut); 
    
    % BUILD TEST DATA (WITH APPROPRIATE FREUQENCIES
    G        = fft(data,[],2);
    G0       = zeros(size(data));
    G0(:,fn) = G(:,fn);
    g0       = 2*real(ifft(G0,[],2));
    g        = g0(:,1:time_cut); 
    
    % BUILD TRAIN DATA (WITH APPROPRIATE FREUQENCIES
    if cortest || scltest
        W        = fft(train_data,[],2);
        W0       = zeros(size(data));
        W0(:,fn) = W(:,fn);
        w0       = 2*real(ifft(W0,[],2));
        w        = w0(:,1:time_cut);     
    end
        
    % PERFORM TEMPERATURE COMPENSATION
    if scltest
        z = zeros(size(g));
        for nn = 1:size(g,1)
            strch = findstretch(g(nn,:).', w(nn,:).');
            z(nn,:) = circscale(g(nn,:).', strch).';
        end
    end
    
    % COMPUTE CORRELATION COEFFICIENTS
    corrs(i,2) = corr(g(:),h(:),'type','Pearson'); ii = 1;
    if cortest, corrs(i,2+ii) = corr(g(:),w(:),'type','Pearson'); ii = ii + 1; end
    if scltest, corrs(i,2+ii) = corr(z(:),w(:),'type','Pearson'); ii = ii + 1; end
    
    fprintf(repmat('\b', 1, 22+8*C));
    fprintf('%04d/%04d(%01d): ',i,size(corrs,1),corrs(i,1)>0);
    fprintf('%6.3f, ', corrs(i,2:end))
end
fprintf('\n');