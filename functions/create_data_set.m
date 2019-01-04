function data_set = create_data_set(num_samplings, num_samples, undamaged_data, damaged_data, train_data, label_undam, label_dam)
%Creates a data_set that can be used for damage detection using dictionary
%learning. The data_set includes a number of samplings, which can each be
%used to create one correlation value for damage detection. Each sampling
%is made up of a certain number of samples from f_test_data_undam and
%f_test_data_dam. The first half of the samplings will be from the
%undamaged data, the second half will be from the damaged data.
%
%Inputs - num_samplings     : Number of samplings to do for each plate pair
%                             (half undamaged, half damaged) 
%         num_samples       : Number of samples to draw from
%                             "f_test_data_undam" or "f_test_data_dam" for
%                             each sampling  
%         f_test_data_undam : file name for the undamaged test data signals
%         f_test_data_dam   : file name for the damaged test data signals
%Output - data_set          : Includes "num_samplings" samplings of
%                             "num_samples" samples each. The plate the the
%                             sampling was drawn from, the label indicating
%                             damage state, the random indices indicating
%                             where the samples were taken from for each
%                             sampling, and the samplings themselves, are
%                             all contained in data_set.

%SET LEARNING PARAMETERS, LOAD DATA FILES
data_set.num_samplings = num_samplings; 
data_set.num_samples = num_samples; 
% undamaged_data = importdata(f_test_data_undam);
% damaged_data = importdata(f_test_data_dam);

%CREATE DATASET
fprintf('Creating dataset... ');
for j = 1:num_samplings
    dim = round(sqrt((num_samples/10000)*length(undamaged_data))); grid_size=[dim,dim];
    random_rows = jitter_sample([100,100], grid_size); 

    if j > ceil(num_samplings/2)
        data_set.(['p',num2str(j)]) = char(label_undam); %what plate was used
        data_set.(['l',num2str(j)]) = 1; %undamaged
        data_set.(['i',num2str(j)]) = random_rows; %what rows indices are used for this sample
        data_set.(['s',num2str(j)]) = undamaged_data(random_rows,:); %actual sample set
        data_set.(['t',num2str(j)]) = train_data(random_rows,:); %actual sample set
    else
        data_set.(['p',num2str(j)]) = char(label_dam); %what plate was used
        data_set.(['l',num2str(j)]) = -1; %damaged
        data_set.(['i',num2str(j)]) = random_rows; %what rows indices are used for this sample
        data_set.(['s',num2str(j)]) = damaged_data(random_rows,:); %actual sample set
        data_set.(['t',num2str(j)]) = train_data(random_rows,:); %actual sample set
    end
end
fprintf('complete\n');