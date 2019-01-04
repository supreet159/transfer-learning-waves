% -------------------------------------------------------------------------
% example_damagedetection.m
% -------------------------------------------------------------------------
% In this script, we demonstrate the use of dictionary learning for damage
% detection. 
%

% -------------------------------------------------------------------------
% Copyright (C) 2016-2017  Joel B. Harley
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or (at 
% your option) any later version. You should have received a copy of the 
% GNU General Public License along with this program. If not, see 
% <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------
% Last updated: November 9, 2017
% -------------------------------------------------------------------------



clear;
addpath('.\functions\');
addpath('..\ksvd\');
addpath('..\Data\Wavefields\');


%% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECIFY PARAMETERS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SET TRAINING DATA AND TEST DATA (NOTE: THIS SCRIPT TRAINS WITH ALL FILES WITH A DIFFERENT MATERIAL / THICKNESS
train_files       = {'AL0625_1_undam', 'AL0625_2_undam', 'ST0625_undam', 'AL25_undam'};         %, 'AL0625_1_dam', 'AL0625_2_dam', 'ST0625_dam', 'AL25_dam'
test_files        = {'AL0625_1_undam', 'AL0625_2_undam', 'ST0625_undam', 'AL25_undam', 'AL0625_1_dam', 'AL0625_2_dam', 'ST0625_dam', 'AL25_dam'};
test_labls        = {'AL-1','AL-2','ST','AL-T'};
excitation_file   = 'currentSignal';                                                            % Files with excitation waveform
test_undam_states = [1 1 2 3 -1 -1 -2 -3];                                                      % States of each undamaged files (negative means damage, values mean material)
mass_loc          = [nan, nan; nan, nan; nan, nan; nan, nan; 36, 67; 23, 80; 22, 82; 34, 79;];  % Mass locations [mm]
mass_radius       = [nan, nan, nan, nan, 8 8 8 8];                                              % Mass radius [mm]
filt_sigma_test   = [2 2 2 4 2 2 2 4];
filt_sigma_train  = [2 2 2 4];

% DICTIONARY LEARNING PARAMETERS (see function learn_complex_dict for details)
fn_train      = 20:1000;                 % Frequency indices to learn from 
DictIter      = 10;                     % Number of dictionary learning iterations
num_atom      = 300;                    % Number of atoms
sparse_learn  = 1;                      % Number of sparse coefficients

% DAMAGE DETECTION PARAMETERS
fn            = 50:1:500;               % Frequency indices to test with
                                        %   10:110 seems to work pretty well
sparse_test   = 1;                      % Number of sparse coefficients
time_cut      = 100;                    % Number of time signals in the rows of the data that should be used for the correlation
Fs            = 300000;                 % Sampling rate

% PLOTTING PARAMETERS
fontsize = 8;                           % Plot fontsize       


%%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% NO NEED TO TRAIN ANYTHING BELOW THIS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% DEFINE CONSTANTS 
Q = 3001;                               % Number of time samples in each measurement
B = 10000;                              % Number of measurements in each file
M = length(test_files);                 % Number of test files
P = length(train_files);                % Number of training files
R = length(fn);                         % Number of frequencies
rows = 1:B;                             % Measurements to use in tests

% EXTRACT EXCITATION SIGNAL
str = load('currentSignal');            % Load excitation
s   = str.currentSignal;                % Assign excitation

% INITIALIZE OUTPUTS
y     = zeros(B,Q,M);                   % Reconstructed data
cc    = zeros(M);                       % Correlation coefficients
ccfrq = zeros(Q,M);                     % Correlation coefficients over freq.
ccmap = zeros(B,M);                     % Correlation coefficients over space

%%
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COLLECT TRAINING DATA
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOOP OVER TRAINING FILES
G = zeros(B,Q,P);
for p = 1:P

    % LOAD TRAINING DATA
    str = load(train_files{p});
    train_data = getfield(str, train_files{p});                  % Get training data
    %train_data(:,time_cut:end) = 0;                             % Cut data in time
    train_data = ifft(fft(train_data.').*conj(fft(s, size(train_data,2)))).';  % Pulse compress
    train_data(:,time_cut:end) = 0;                              % Cut data in time again
    G(:,:,p) = fft(train_data, [], 2);                           % Fourier transform
    
    % FILTER DATA
    Gi = permute(reshape(G(:,:,p).',Q,sqrt(B), sqrt(B)),[2 3 1]);
    Ga = imgaussfilt(real(Gi), filt_sigma_train(p))+1j*imgaussfilt(imag(Gi), filt_sigma_train(p));
    G(:,:,p) = reshape(permute(Ga, [3 1 2]), Q, B).'; 
    
    % NORMALIZE EACH MEASUREMENT TO HAVE UNIT ENERGY
    G(:,:,p) = G(:,:,p)./sqrt(sum(abs(G(:,:,p)).^2,2))*sqrt(Q);
    
end

% CREATE TRAINING DATA SET WITH ONLY SPECIFIED TEST FREQUENCIES
G0 = zeros(B,Q,P);
G0(:,fn,:) = G(:,fn,:); 
g = 2*real(ifft(G0, [], 2));


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COLLECT TEST DATA
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOOP OVER TEST FILES
H = zeros(B,Q,M);
for m = 1:M

    % GET TEST DATA 
    str = load(test_files{m});                                
    f_test_data_undam = getfield(str, test_files{m});           % Get training data
    %f_test_data_undam(:,time_cut:end) = 0;                      % Cut data in time
    f_test_data_undam = ifft(fft(f_test_data_undam.').*conj(fft(s, size(f_test_data_undam,2)))).';  % Pulse compress
    f_test_data_undam(:,time_cut:end) = 0;                      % Cut data in time
    H(:,:,m) = fft(f_test_data_undam, [], 2); 
    
    % FILTER DATA
    Hi = permute(reshape(H(:,:,m).',Q,sqrt(B), sqrt(B)),[2 3 1]);
    Ha = imgaussfilt(real(Hi), filt_sigma_test(m))+1j*imgaussfilt(imag(Hi), filt_sigma_test(m));
    H(:,:,m) = reshape(permute(Ha, [3 1 2]), Q, B).'; 
    
    % NORMALIZE EACH MEASUREMENT TO HAVE UNIT ENERGY
    H(:,:,m) = H(:,:,m)./sqrt(sum(abs(H(:,:,m)).^2,2))*sqrt(Q);

end

% CREATE TEST DATA SET WITH ONLY SPECIFIED FREQUENCIES
H0 = zeros(B,Q,M);
H0(:,fn,:) = H(:,fn,:); 
h = 2*real(ifft(H0, [], 2));


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RECONSTRUCT DATA
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOOP OVER TEST FILES
for m = 1:M

    % DETERMINE FILES TO TRAIN WITH
    testf = strsplit(test_files{m}, '_', 'CollapseDelimiters', true);
    testf = testf{1};
    tranf = cellfun(@(x) strsplit(x, '_', 'CollapseDelimiters', true), train_files, 'UniformOutput', false);
    tranf = tranf(:); tranf = cellfun(@(x) x{1}, tranf, 'UniformOutput', false);
    usefiles = cellfun(@(x) ~strcmp(x, testf), tranf);
    
    % TRAIN DICTIONARY
    Di = randn(B,num_atom)+1j*randn(B,num_atom);                                  % Initial Dictionary (not normalized)
    [ D, V ] = ksvd_simple( ...                                                   % Learned dictionary
                reshape(G(:,fn_train,usefiles), B, length(fn_train)*sum(usefiles)), ...
                Di, sparse_learn, DictIter ...
               );  

    % RECONSTRUCT WAVEFIELD
    [~,yu] = reconstruct(h(rows,:,m).',D,rows,fn,sparse_test); 
    y(:,:,m) = yu.';
           
    % COMPUTE CORRELATION COEFFICIENTS
    cc(m)      = sum(sum(y(:,:,m).*h(:,:,m)))/norm(y(:,:,m),'fro')/norm(h(:,:,m),'fro');    
    ccfrq(:,m) = sum(fft(y(:,:,m),[],2).*conj(fft(h(:,:,m),[],2)),1)./sqrt(sum(abs(fft(y(:,:,m),[],2)).^2,1))./sqrt(sum(abs(fft(h(:,:,m),[],2)).^2,1));
    ccmap(:,m) = sum(y(:,:,m).*h(:,:,m),2)./sqrt(sum(abs(y(:,:,m)).^2,2))./sqrt(sum(abs(h(:,:,m)).^2,2));

end


%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RECONSTRUCT DATA
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SET PARAMETERS
ttm = [26 28 32 20 26 28 31 20];        % Times to plot wavefields
m0  = [3 3 1 3 3 3 1 3];                % Surrogate measurements
limmin = -0.2;                          % Lower bound for data amplitude
limmax = 0.2;                           % Upper bound for data amplitude

% SET DATA
m1 = 1;                                 % Surrogate measurements
m2 = 7;                                 % Test measurment
ii = ttm(m2);                           % Test time
sur_data = reshape(h(:,ii,m1),100,100); % Surrogate data
dam_data = reshape(h(:,ii,m2),100,100); % Test data

% PLOT WAVEFIELDS
figure(124)
set(gcf, 'Units', 'Inches', 'Position', [1 1 6.5 1.8 ])

ax(1) = subplot(131); 
imagesc(sur_data, [limmin limmax]); 
set(gca, 'fontsize', fontsize)
xlabel({'Length (mm)', '(a) Surrogate Plate'}, 'fontsize', fontsize)
ylabel('Width (mm)', 'fontsize', fontsize)
axis square;
axis xy;
colormap('gray')

ax(2) = subplot(132); 
imagesc(dam_data, [limmin limmax]); 
set(gca, 'fontsize', fontsize)
xlabel({'Length (mm)', '(b) Test Plate (Mass)'}, 'fontsize', fontsize)
ylabel('Width (mm)', 'fontsize', fontsize)
axis square;
hold on;
circle2(mass_loc(m2,1), mass_loc(m2,2), mass_radius(m2))
hold off;
axis xy;
colormap('gray')

ax(3) = subplot(133); 
imagesc(dam_data-sur_data, [limmin limmax]); 
set(gca, 'fontsize', fontsize)
xlabel({'Length (mm)', '(c) Subtraction'}, 'fontsize', fontsize)
ylabel('Width (mm)', 'fontsize', fontsize)
axis square;
axis xy;
colormap('gray')


hbar=colorbar;
set(hbar, 'Position', [.88 .31 .0281 .65])
for i=1:3
      pos=get(ax(i), 'Position');
      set(ax(i), 'Position', [pos(1)-0.15*pos(3) pos(2)+0.15*pos(4) 0.825*pos(3) pos(4)]);
end
drawnow;
    
print('surrogate_example.eps', '-depsc')


%%

% SET PARAMTERS
tme0 = 30;          % Time to plot

% PLOT WAVEFIELDS
figure(10)
set(gcf, 'Units', 'Inches', 'Position', [1 1 3.5 3.75 ])
for m = 1:4    

    mydata = reshape(h(:,tme0,m),100,100);
    
    subplot(2,2,m)
    imagesc(mydata, [limmin limmax])
    set(gca, 'fontsize', fontsize)
    xlabel({'Length (mm)', ['(' char(96+m) ') ' test_labls{m}]}, 'fontsize', fontsize)
    ylabel('Width (mm)', 'fontsize', fontsize)
    axis square;
    axis xy;
    colormap('gray')
    
end
print('./figures/all_wavefield_undamaged.eps', '-depsc')
print('./figures/all_wavefield_undamaged.png', '-dpng')

figure(11)
set(gcf, 'Units', 'Inches', 'Position', [1 1 3.5 3.75 ])
for m = 1:4    

    mydata = reshape(h(:,tme0,m+4),100,100);
    
    subplot(2,2,m)
    imagesc(mydata, [limmin limmax])
    set(gca, 'fontsize', fontsize)
    xlabel({'Length (mm)', ['(' char(96+m) ') ' test_labls{m}]}, 'fontsize', fontsize)
    ylabel('Width (mm)', 'fontsize', fontsize)
    axis square;
    axis xy;
    hold on;
    circle2(mass_loc(m+4,1), mass_loc(m+4,2), mass_radius(m+4))
    hold off;
    colormap('gray')
    
end
print('./figures/all_wavefield_damaged.eps', '-depsc')
print('./figures/all_wavefield_damaged.png', '-dpng')

figure(12)
set(gcf, 'Units', 'Inches', 'Position', [1 1 3.5 3.75 ])
for m = 1:4    

    mydata = reshape(y(:,tme0,m+4),100,100);
    
    subplot(2,2,m)
    imagesc(mydata, [limmin limmax])
    set(gca, 'fontsize', fontsize)
    xlabel({'Length (mm)', ['(' char(96+m) ') ' test_labls{m}]}, 'fontsize', fontsize)
    ylabel('Width (mm)', 'fontsize', fontsize)
    axis square;
    axis xy;
    %hold on;
    %circle2(mass_loc(m+4,1), mass_loc(m+4,2), mass_radius(m+4))
    %hold off;
    colormap('gray')
    
end
print('./figures/all_wavefield_reconstructed.eps', '-depsc')
print('./figures/all_wavefield_reconstructed.png', '-dpng')

figure(13)
set(gcf, 'Units', 'Inches', 'Position', [1 1 3.5 3.75 ])
for m = 1:4    

    mydata = reshape(h(:,tme0,m)-y(:,30,m),100,100);
    
    subplot(2,2,m)
    imagesc(mydata, [limmin limmax])
    set(gca, 'fontsize', fontsize)
    xlabel({'Length (mm)', ['(' char(96+m) ') ' test_labls{m}]}, 'fontsize', fontsize)
    ylabel('Width (mm)', 'fontsize', fontsize)
    axis square;
    axis xy;
    %hold on;
    %circle2(mass_loc(m+4,1), mass_loc(m+4,2), mass_radius(m+4))
    %hold off;
    colormap('gray')
    
end
print('./figures/all_wavefield_subtracted_undamaged.eps', '-depsc')
print('./figures/all_wavefield_subtracted_undamaged.png', '-dpng')

figure(14)
set(gcf, 'Units', 'Inches', 'Position', [1 1 3.5 3.75 ])
for m = 1:4    

    mydata = reshape(h(:,tme0,m+4)-y(:,30,m+4),100,100);
    
    subplot(2,2,m)
    imagesc(mydata, [limmin limmax])
    set(gca, 'fontsize', fontsize)
    xlabel({'Length (mm)', ['(' char(96+m) ') ' test_labls{m}]}, 'fontsize', fontsize)
    ylabel('Width (mm)', 'fontsize', fontsize)
    axis square;
    axis xy;
    hold on;
    circle2(mass_loc(m+4,1), mass_loc(m+4,2), mass_radius(m+4))
    hold off;
    colormap('gray')
    
end
print('./figures/all_wavefield_subtracted_damaged.eps', '-depsc')
print('./figures/all_wavefield_subtracted_damaged.png', '-dpng')

figure(15)
set(gcf, 'Units', 'Inches', 'Position', [1 1 3.5 3.75 ])
for m = 1:4    

    mydata = reshape(h(:,tme0,m+4)-g(:,30,m0(m)),100,100);
    
    subplot(2,2,m)
    imagesc(mydata, [limmin limmax])
    set(gca, 'fontsize', fontsize)
    xlabel({'Length (mm)', ['(' char(96+m) ') ' test_labls{m}]}, 'fontsize', fontsize)
    ylabel('Width (mm)', 'fontsize', fontsize)
    axis square;
    axis xy;
    hold on;
    circle2(mass_loc(m+4,1), mass_loc(m+4,2), mass_radius(m+4))
    hold off;
    colormap('gray')
    
end
print('./figures/all_wavefield_surrsubtracted_damaged.eps', '-depsc')
print('./figures/all_wavefield_surrsubtracted_damaged.png', '-dpng')

%%

% SET PARAMTERS
spcindx = sub2ind([sqrt(B) sqrt(B)], 67, 39);       % Plot location
t = 1/Fs:1/Fs:Q/Fs;                                 % Time axis

% PLOT TIME SIGNALS
figure(20)
set(gcf, 'Units', 'Inches', 'Position', [1 1 3.5 2.4 ], 'Name', 'Time Data (No Damage)', 'color', [1 1 1])
set(gca, 'fontsize', fontsize)
for m = 1:4    

    mydata = h(spcindx,:,m);
    
    subplot(2,2,m);
    plotpretty(@plot, t*1000, mydata, 'fontsize', fontsize, 'axis', [0 t(time_cut)*1.5*1000 limmin limmax]);
    axis([0 t(time_cut)*1.5*1000 limmin limmax])
    xlabel({'Time (ms)', ['(' char(96+m) ') ' test_labls{m}]}, 'fontsize', fontsize)
    if mod(m,2), ylabel('Amplitude', 'fontsize', fontsize); end
    
end
export_fig('./figures/time_undamaged.eps', '-eps', '-png')

figure(21)
set(gcf, 'Units', 'Inches', 'Position', [1 1 3.5 2.4 ], 'Name', 'Time Data (Damage)', 'color', [1 1 1])
for m = 1:4    

    mydata = h(spcindx,:,m+4);
    
    subplot(2,2,m)
    plotpretty(@plot, t*1000, mydata, 'fontsize', fontsize, 'axis', [0 t(time_cut)*1.5*1000 limmin limmax]);
    axis([0 t(time_cut)*1.5*1000 limmin limmax])
    xlabel({'Time (ms)', ['(' char(96+m) ') ' test_labls{m}]}, 'fontsize', fontsize)
    if mod(m,2), ylabel('Amplitude', 'fontsize', fontsize); end
    
end
export_fig('./figures/time_damaged.eps', '-eps', '-png')

figure(22)
set(gcf, 'Units', 'Inches', 'Position', [1 1 3.5 2.4 ], 'Name', 'Recon. Time Data (Damage)', 'color', [1 1 1])
for m = 1:4    

    mydata = y(spcindx,:,m+4);
    
    subplot(2,2,m)
    plotpretty(@plot, t*1000, mydata, 'fontsize', fontsize, 'axis', [0 t(time_cut)*1.5*1000 limmin limmax]);
    axis([0 t(time_cut)*1.5*1000 limmin limmax])
    xlabel({'Time (ms)', ['(' char(96+m) ') ' test_labls{m}]}, 'fontsize', fontsize)
    if mod(m,2), ylabel('Amplitude', 'fontsize', fontsize); end
    
end
export_fig('./figures/time_recon_damage.eps', '-eps', '-png')

figure(23)
set(gcf, 'Units', 'Inches', 'Position', [1 1 3.5 2.4 ], 'Name', 'Subtracted Time Data (No Damage)', 'color', [1 1 1])
for m = 1:4    

    mydata = y(spcindx,:,m) - h(spcindx,:,m);
    
    subplot(2,2,m)
    plotpretty(@plot, t*1000, mydata, 'fontsize', fontsize, 'axis', [0 t(time_cut)*1.5*1000 limmin limmax]);
    axis([0 t(time_cut)*1.5*1000 limmin limmax])
    xlabel({'Time (ms)', ['(' char(96+m) ') ' test_labls{m}]}, 'fontsize', fontsize)
    if mod(m,2), ylabel('Amplitude', 'fontsize', fontsize); end
    
end
export_fig('./figures/time_subtracted_undamaged.eps', '-eps', '-png')

figure(24)
set(gcf, 'Units', 'Inches', 'Position', [1 1 3.5 2.4 ], 'Name', 'Subtracted Time Data (Damage)', 'color', [1 1 1])
for m = 1:4    

    mydata = y(spcindx,:,m+4) - h(spcindx,:,m+4);
    
    subplot(2,2,m)
    plotpretty(@plot, t*1000, mydata, 'fontsize', fontsize, 'axis', [0 t(time_cut)*1.5*1000 limmin limmax]);
    xlabel({'Time (ms)', ['(' char(96+m) ') ' test_labls{m}]}, 'fontsize', fontsize)
    if mod(m,2), ylabel('Amplitude', 'fontsize', fontsize); end
    
end
export_fig('./figures/time_subtracted_damaged.eps', '-eps', '-png')
%%

% SET PARAMTERS
m = 2;                  % Measurement to plot

% PLOT TIME SIGNALS
figure(30)
set(gcf, 'Units', 'Inches', 'Position', [1 1 3.5 4 ], 'Name', 'Time Data (No Damage)', 'color', [1 1 1])
set(gca, 'fontsize', fontsize)    

mydata1 = h(spcindx,:,m);
mydata2 = y(spcindx,:,m);
mydata3 = g(spcindx,:,m0(m));
mydata4 = h(spcindx,:,m+4);
mydata5 = y(spcindx,:,m+4);

subplot(3,1,1);
h0 = plotpretty(@plot, t*1000, [mydata3(:) mydata1(:)], 'fontsize', fontsize, 'axis', [0 t(time_cut)*1.5*1000 limmin limmax], 'colorset', 'ColorBrewerBluesTwo');
set(h0(2), 'linewidth', 2)
set(h0(1), 'linewidth', 0.5)
axis([0 t(time_cut)*2*1000 limmin-0.025 limmax+0.025])
xlabel({'Time (ms)', ['(' char(96+1) ') AL-2 Test (No Mass) / ST Surrogate' ]}, 'fontsize', fontsize)
ylabel('Amplitude', 'fontsize', fontsize);
line([0.38 0.41], [0.18 0.18], 'color', [166,189,219]/255, 'linewidth', 2')
line([0.38 0.41], [0.08 0.08], 'color', [4  ,90 ,141]/255, 'linewidth', 1')
ht(1) = text(0.42 , 0.18, 'ST Surrogate', 'fontsize', fontsize-1, 'color', [166,189,219]/255/2);
ht(2) = text(0.42 , 0.08, 'AL-2 Test (No Mass)', 'fontsize', fontsize-1, 'color', [4  ,90 ,141]/255/2);


subplot(3,1,2);
h0 = plotpretty(@plot, t*1000, [mydata2(:) mydata1(:)], 'fontsize', fontsize, 'axis', [0 t(time_cut)*1.5*1000 limmin limmax], 'colorset', 'ColorBrewerBluesTwo');
set(h0(2), 'linewidth', 2)
set(h0(1), 'linewidth', 0.5)
axis([0 t(time_cut)*2*1000 limmin-0.025 limmax+0.025])
xlabel({'Time (ms)', ['(' char(96+2) ') AL-2 Test (No Mass) / AL-2 Reconstruction']}, 'fontsize', fontsize)
ylabel('Amplitude', 'fontsize', fontsize); 
line([0.38 0.41], [0.18 0.18], 'color', [166,189,219]/255, 'linewidth', 2')
line([0.38 0.41], [0.08 0.08], 'color', [4  ,90 ,141]/255, 'linewidth', 1')
ht(3) = text(0.42 , 0.18, 'AL-2 Reconstruction', 'fontsize', fontsize-1, 'color', [166,189,219]/255/2);
ht(4) = text(0.42 , 0.08, 'AL-2 Test (No Mass)', 'fontsize', fontsize-1, 'color', [4  ,90 ,141]/255/2);

subplot(3,1,3);
h0 = plotpretty(@plot, t*1000, [mydata5(:) mydata4(:)], 'fontsize', fontsize, 'axis', [0 t(time_cut)*1.5*1000 limmin limmax], 'colorset', 'ColorBrewerBluesTwo');
set(h0(2), 'linewidth', 2)
set(h0(1), 'linewidth', 0.5)
axis([0 t(time_cut)*2*1000 limmin-0.025 limmax+0.025])
xlabel({'Time (ms)', ['(' char(96+3) ') AL-2 Test (With Mass) / AL-2 Reconstruction']}, 'fontsize', fontsize)
ylabel('Amplitude', 'fontsize', fontsize); 
line([0.38 0.41], [0.18 0.18], 'color', [166,189,219]/255, 'linewidth', 2')
line([0.38 0.41], [0.08 0.08], 'color', [4  ,90 ,141]/255, 'linewidth', 1')
ht(5) = text(0.42 , 0.18, 'AL-2 Reconstruction', 'fontsize', fontsize-1, 'color', [166,189,219]/255/2, 'fontname', 'Calibri');
ht(6) = text(0.42 , 0.08, 'AL-2 Test (With Mass)', 'fontsize', fontsize-1, 'color', [4  ,90 ,141]/255/2, 'fontname', 'Calibri');

for ii = 1:length(ht)
    set(ht(ii), 'fontsize', fontsize-1)
end



export_fig('./figures/example_surr_recon.eps', '-eps', '-png')
