function [accuracy_kmeans, thresh_kmeans, margin] = classify_correlations(corrs, if_plot, if_print,num_samplings)
%Given the correlations produced using the "recovery_correlation" function
%(including corresponding labels of damage state), an optimal threshold and
%a k-means threshold are produced. The accuracy, threshold, and margin are
%printed and a plot of the undamaged and damaged correlation values, and
%the thresholds, is produced. 
%
%Inputs -  corrs     : Correlation outcomes from the "recovery_correlation"
%                      function that contains labels and correlations 
%          if_plot   : If true, then the correlations will be plotted
%          if_print  : If true, then the results will be printed
%Outputs - accuracy  : The best possible accuracy for the given
%                      correlations and labels 
%          threshold : Median threshold that gives the lowest error
%          margin    : The distance between the lowest undamaged
%                      correlationa and the highest damaged correlation
%                      (full scale is 0% to 100%) 


% SET RANDOM NUMBER SEED
rng(1, 'twister') 

% SPECIFY CONSTANTS
M = size(corrs,1);
C = size(corrs,2)-1;
fontsize = 8;

% INITIALIZE OUTPUT
thresh_kmeans = zeros(C,1);
accuracy_kmeans = zeros(C,1);
margin = zeros(C,1);

% LOOP OVER CORRELATION COEFFICIENTS
for ii = 1:C
    
    % COMPUTE K-MEANS THRESHOLD 
    [~, C0] = kmeans(corrs(:,1+ii), 2);
    thresh_kmeans(ii) = (C0(1)+C0(2))/2;

    % COMPUTE K-MEANS ACCURACY
    accuracy_kmeans(ii) = mean([corrs(corrs(:,1+ii) > thresh_kmeans(ii),1)== 1; ...
                                corrs(corrs(:,1+ii) < thresh_kmeans(ii),1)==-1]);
                        
    % COMPUTE MARGIN
    min_corr   = min(corrs((corrs(:,1)== 1), 1+ii));
    max_corr   = max(corrs((corrs(:,1)==-1), 1+ii));
    margin(ii) = max(min_corr - max_corr   , 0);


    %PLOT RESULTS
    if if_plot
        figure(123);  
        set(gcf, 'Units', 'Inches', 'Position', [1 1 8 3])
        subplot(1, C, ii)
        plot([0,M/2],[thresh_kmeans(ii),thresh_kmeans(ii)],'--k','linewidth',1);
        hold on;
        plot(corrs(corrs(:,1)==1,1+ii),'o', 'MarkerFaceColor', [158,202,225]/255, 'color', [1 1 1]);
        plot(corrs(corrs(:,1)==-1,1+ii),'x', 'MarkerFaceColor', [49,130,189]/255, 'color', [49,130,189]/255);
        %axis([0 M/2 0.3 1])
        title(['Accuracy: ',num2str(100*accuracy_kmeans(ii)),'%, ','Margin: ',num2str(100*margin(ii)),'%']); 
        ylabel('Recovery Correlation');
        xlabel('Monte Carlo Trials');
        hh = legend('Kmeans Threshold','Undamaged','Damaged','Location','southoutside','Orientation', 'horizontal');
        legend('boxoff')
        set(hh, 'fontsize', fontsize)
        set(gca,'fontsize', fontsize)
        hold off
        drawnow;
    end
    
    %PRINT RESULTS
    if if_print
        fprintf('Kmeans Threshold -  Accuracy: %3.0f%%, Threshold: %3.0f%% \n',accuracy_kmeans(ii)*100,thresh_kmeans(ii)*100);
    end


end

fprintf('\n')

end

