%%% KMeans function

function [KMeans_res]= kmeans_fun(esl,num_clusters,Replicates,lonlat)

% esl= extreme surge levels
% num_clusters= number of clusters
% Replicates= for kmeans

%% Evaluate the clustering using HandMade Criterion

rng(1);
opts = statset('Display','off');

[idx,C,~,D]= kmeans(esl,num_clusters,'Distance','correlation',...
    'Replicates',Replicates,'Options',opts,'OnLinePhase','on');

ref_gauges= nan(num_clusters,3);

ref_series_position= nan(num_clusters,1);

% reference series
for i= 1: num_clusters
    
    % centroid
    fis= find(abs(D(:,i))== min(abs(D(:,i))));
    fis= fis(1);
    ref_gauges(i,:)= [lonlat(fis,1),lonlat(fis,2),idx(fis)];
    ref_series_position(i)= fis;
    
end

% sort by latitude
[~,ref_series_sort]= sort(ref_gauges(:,2));
ref_gauges= ref_gauges(ref_series_sort,:);

coordinates_idx= [lonlat, nan(length(lonlat),1)];

for i= 1: num_clusters
   
    coordinates_idx(idx== ref_series_sort(i),3)= ref_series_sort(i);
  
end


% Save the idx of this combination
KMeans_res{2,1}= coordinates_idx;
KMeans_res{2,2}= num_clusters;
KMeans_res{2,3}= D;
KMeans_res{2,4}= C;
KMeans_res{2,5}= ref_gauges;
KMeans_res{2,6}= ref_series_position;

KMeans_res{1,1}= 'gauge indexs idx';
KMeans_res{1,2}= 'number of clusters used';
KMeans_res{1,3}= 'distance to the centers';
KMeans_res{1,4}= 'theoretical centroids';
KMeans_res{1,5}= 'reference gauges and idx';
KMeans_res{1,6}= 'reference series positions';



