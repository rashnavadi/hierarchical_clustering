Here is why we see only 30 observation indexes:
https://www.mathworks.com/help/stats/dendrogram.html
Dendrogram plot - MATLAB dendrogram
This MATLAB function generates a dendrogram plot of the hierarchical binary cluster tree.
www.mathworks.com

dendrogram(tree) generates a dendrogram plot of the hierarchical binary cluster tree. A dendrogram consists of many U-shaped lines that connect data points in a hierarchical tree. The height of each U represents the distance between the two data points being connected.
If there are 30 or fewer data points in the original data set, then each leaf in the dendrogram corresponds to one data point.
If there are more than 30 data points, then dendrogram collapses lower branches so that there are 30 leaf nodes. As a result, some leaves in the plot correspond to more than one data point.

And below is the code that i used to run hierarchical clustering to get those outputs, however, before running the code, i needed to prepare the Mtx matrix myself by concatenating the FC matrices, for example, for left FLE with six subjects, each having 140 timepoints, for DMN for example: for each subject i have 140 rows of (38x38 -38) datapoints/columns, concatenating these vectors (i flattened the FC matrices and turned them to vectors) for all six subjects i obtained: 6 (subjects) x 140 (timepoints) rows/observations ( = 840) and (38x38 -38) or 1406 columns or features for my input data to hierarchical clustering, since i have 840 observations so matlab cut it down to 30 as explained above in their website.
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
