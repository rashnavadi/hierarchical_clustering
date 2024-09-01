%% hierarchical clustering on Left FLE data that seems to have only one state.
% written by Tahereh Rashnavadi 
% june 2024
close all 
clear
% Number of subjects and time points
% num_subjects = 5; % left fle, leave one out (LOO)
num_subjects = 6; % left fle
% num_subjects = 9; % healthy control
% num_subjects = 10; % right fle
% num_subjects = 25; % all fle
num_timepoints = 140;

% Total number of observations
total_observations = num_subjects * num_timepoints;

% Create an array to hold subject labels
subject_labels = [];

% Populate the subject_labels array
for subj = 1:num_subjects
    subject_labels = [subject_labels; repmat(subj, num_timepoints, 1)];
end

cd /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/FLE_2024/all_subjs_for_paper/hierarchical_clustering

% Now load your matrix Mtx
% Assuming Mtx is loaded here with size [total_observations, DMN: 38x38]
% Assuming Mtx is loaded here with size [total_observations, SN: 12x12]




% the output from HOCo analysis
% left FLE, sensorimotor
% load /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/FLE_2024/all_subjs_for_paper/left_FLE/kmeans_output/sensorimotor/reordered_ROIs/k_means_output_2.mat
% LOO left FLE, sensorimotor
% load /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/FLE_2024/all_subjs_for_paper/left_FLE/leave_one_out/leave_one_out_kmeans/sensori/reordred/k_means_output_2.mat
% left FLE, DMN
load /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/FLE_2024/all_subjs_for_paper/left_FLE/kmeans_output/DMN/reordered_ROIs/k_means_output_2.mat
% LOO left FLE, DMN
% load /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/FLE_2024/all_subjs_for_paper/left_FLE/leave_one_out/leave_one_out_kmeans/DMN/reordered/k_means_output_2.mat

% HC, sensorimotor
% load /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/FLE_2024/all_subjs_for_paper/controls/kmeans_output/sensorimotor/reordered_ROIs/reordered_k_means_output_2.mat
% HC, DMN
% load /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/FLE_2024/all_subjs_for_paper/controls/kmeans_output/DMN/my_reordering_leftright/k_means_output_2.mat

% right fle, sensorimotor
% load /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/FLE_2024/all_subjs_for_paper/right_FLE/kmeans_output/sensorimotor/reordered_ROIs/reordered_k_means_output_2.mat
% right fle, DMN
% load /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/FLE_2024/all_subjs_for_paper/right_FLE/kmeans_output/DMN/my_reordering_leftright/k_means_output_2.mat

% all fle, sensorimotor
% load /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/FLE_2024/all_subjs_for_paper/all_subjects/kmeans_output/sensorimotor/reordered_ROIs/k_means_output_2.mat
% all fle, DMN
% load /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/FLE_2024/all_subjs_for_paper/all_subjects/kmeans_output/DMN/reordered_ROIs/k_means_output_2.mat

% Compute the distance matrix using the cityblock metric
distances = pdist(Mtx, 'cityblock');

% Perform hierarchical clustering
Z = linkage(distances, 'ward');


% Cut the dendrogram to form clusters
num_clusters = 2; % Specify the number of clusters
clusters = cluster(Z, 'maxclust', num_clusters);

% Determine the dominant subject in each cluster
dominant_subject = zeros(num_clusters, 1);
for i = 1:num_clusters
    % Find the subject labels for all observations in this cluster
    cluster_subjects = subject_labels(clusters == i);
    
    % Determine which subject is most frequent in this cluster
    dominant_subject(i) = mode(cluster_subjects);
end

% Display the dendrogram with cluster labels based on the dominant subject
figure;
H = dendrogram(Z, 0);

% Color the dendrogram branches according to the dominant subject
for i = 1:length(H)
    cluster_index = clusters(str2double(get(H(i), 'UserData')));
    color_index = dominant_subject(cluster_index);
    set(H(i), 'Color', lines(num_clusters)(color_index, :));
end

set(H, 'LineWidth', 2);
title('Hierarchical Clustering Dendrogram');
xlabel('Observation Index');
ylabel('Distance');

% Annotate the dendrogram with subject labels
xticklabels(arrayfun(@(x) sprintf('Subject %d', dominant_subject(x)), clusters, 'UniformOutput', false));

% Optional: Display clusters and their dominant subject
disp('Cluster assignments and dominant subjects:');
for i = 1:num_clusters
    fprintf('Cluster %d: Dominant Subject = %d\n', i, dominant_subject(i));
end


%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%% Assuming 'clusters' is the output from hierarchical clustering and is of size 840x1
% Plot clusters for each subject in separate figures
for subj = 1:num_subjects
    % Extract the clusters for the current subject
    subject_clusters = clusters((subj-1)*num_timepoints + 1 : subj*num_timepoints);
    
    % Create a new figure for each subject
    figure;
    plot(1:num_timepoints, subject_clusters, 'LineWidth', 4);
    yticks([1,num_clusters]);
    ylim([0.5, num_clusters + .5]);
    title(sprintf('Subject %d', subj));
    xlabel('Time Points', 'FontSize', 14);
    ylabel('State Number', 'FontSize', 14);
end
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

%% cut the dendogram based on the number of clusters that you want not the distance between the data
% Cut the dendrogram to form the specified number of clusters
% Find the maximum distance within each of the clusters
% cut_height = max(Z(find(diff(clusters) ~= 0), 3)); % 3 stands for the third column of Z
% Generate the full dendrogram and highlight the specified number of clusters
% Determine the height threshold for the desired number of clusters
% Sort the heights in ascending order
sortedHeights = sort(Z(:,3));

% The cut height is the one that gives us the desired number of clusters
cut_height = sortedHeights(end - num_clusters + 2);
figure;
[H, T, outperm] = dendrogram(Z, 0, 'ColorThreshold', cut_height);
set(H, 'LineWidth', 2); % Adjust the line width after plotting

title('Hierarchical Clustering Dendrogram');
xlabel('Observation Index', 'FontSize', 14);
ylabel('Distance', 'FontSize', 14);

% Highlight the clusters with different colors
colors = lines(num_clusters);
for i = 1:num_clusters
    set(H(T == i), 'Color', colors(i,:));
end
%% cut the dendogram based on the distance between the data

% Define the cut height (threshold) for the dendrogram
% cut_height = 280; % for sensorimotor Adjust this value as needed
% cut_height = 250; % LOO, for sensorimotor Adjust this value as needed
cut_height = 4300; % LOO, for DMN Adjust this value as needed

% cut_height = 800; % for sensorimotor all fles
% cut_height = 9000; % for sensorimotor all fles

% cut_height = 5500; % for DMN, left/right FLE
% cut_height = 8100; % for DMN, HC

% 
% % Plot the dendrogram with a cut at the specified height
figure;
[H, T] = dendrogram(Z, 'ColorThreshold', cut_height);
set(H, 'LineWidth', 2); % Adjust the line width after plotting
title('Hierarchical Clustering Dendrogram (Cut)', 'FontSize', 14);
xlabel('Observation Index', 'FontSize', 14);
ylabel('Distance', 'FontSize', 14);

% Optional: Highlight the cut
line([0, length(Z) + 1], [cut_height, cut_height], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

%% plot the number of clusters vs. threshold of distance in dendrogram
% Initialize variables to store results
maxClusters = num_subjects; % Set maximum number of clusters to 6
numSteps = size(Z, 1); % Number of steps in hierarchical clustering
last_Steps = numSteps - maxClusters -1 + 2:numSteps; % Indices of the last steps
mergeDistances = Z(last_Steps, 3); % Extract the last 6 merge distances

% Number of clusters (from 6 to 1 as we move up the dendrogram)
numClusters = 2:maxClusters;

mergeDistances = flip(mergeDistances);
% Plot the number of clusters vs. merge distances
figure;
plot(numClusters, mergeDistances(1:end-1), '-o', 'LineWidth', 2);
xticks(numClusters); % Set x-axis ticks to be 1 to 6 only

xlabel('Number of Clusters');
ylabel('Distance where Branch is Added/Divided');
title('Distance vs. Number of Clusters');
grid on;


