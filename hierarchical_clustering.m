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

cd /Users/trashnavadi/Documents/Data_Analysis/2023/analyses/FLE_2024/all_subjs_for_paper/hierarchical_clustering
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

% Plot the dendrogram
figure;
H = dendrogram(Z, 0); % Plot the dendrogram and capture handles
% Adjust the line width of the dendrogram
set(H, 'LineWidth', 2); % Adjust the line width as neededtitle('Hierarchical Clustering Dendrogram');
xlabel('Observation Index');
ylabel('Distance');

% Optional: Cut the dendrogram to form clusters
num_clusters = 2; % Specify the number of clusters
clusters = cluster(Z, 'maxclust', num_clusters);

% Display the clusters
disp('Clusters:');
% disp(clusters);

% Perform PCA
[coeff, score] = pca(Mtx);

% Plot the first three principal components in a 3D scatter plot
hFig = figure();
axh = axes('Parent', hFig);
hold(axh, 'all');
unique_clusters = unique(clusters);

% Define a custom colormap with distinct colors for each cluster
custom_colors = lines(num_clusters); % You can use other colormaps or define your own colors

% Plot each cluster with a different color
for i = 1:numel(unique_clusters)
    % Filter data points belonging to cluster i
    cluster_indices = clusters == unique_clusters(i);
%     color = rand(1,3);  % Random color for each cluster
    color = custom_colors(i, :);
    scatter3(axh, score(cluster_indices,1), score(cluster_indices,2), score(cluster_indices,3), 50, color, 'filled');
    hold(axh, 'on');
end

% grid(axh, 'on') % Uncomment this line if you want the grid to be displayed
view(axh, 3); % Ensure the view is in 3D

% Create legend for each cluster
legend_labels = cellstr(num2str((1:num_clusters)', 'Cluster %d'));
legend(axh, legend_labels, 'FontSize', 12); % Adjust the font size as needed

title(axh, '3D PCA Scatter Plot with Clusters', 'FontSize', 14);
xlabel(axh, 'First Principal Component', 'FontSize', 14);
ylabel(axh, 'Second Principal Component', 'FontSize', 14);
zlabel(axh, 'Third Principal Component', 'FontSize', 14);
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
% Assuming 'clusters' is the output from hierarchical clustering and is of size 840x1
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


