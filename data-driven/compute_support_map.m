function [s,dist] = compute_support_map(idx, meanV, F, t, dist_pre, smooth_min_dist, smooth_max_dist)

h = mean(mean(edge_lengths(meanV,F))); % average edge length

%dist = heat_geodesic(meanV,F,idx,t,'Precomputation',dist_pre,'Legacy',false);
% this is bug in Alec's implementation.

%dist = pdist2(meanV,meanV(idx,:),'euclidean');

%[dist,~] = heat_geodesic_naive(meanV,F,idx);
[dist,~] = heat_geodesic_keenan(meanV,F,idx);%,'non_right_mass');

if( max(dist)-min(dist)<0.01*h || max(dist)==min(dist) )
   warning(['Mesh might not be fully connected at vertex',num2str(idx),',used pdist2 instead!\n']);
   dist = pdist2(meanV,meanV(idx,:),'euclidean');
end

td = dist;

td = max(td,smooth_min_dist);
td = min(td,smooth_max_dist);

s = (td - smooth_min_dist) / (smooth_max_dist - smooth_min_dist);