function [output gmm IM]= infer_graph_from_datastream_mmVB(DATA,prior_on_K)

agents = DATA(:,2);

IM = index_manager(agents);
DATA(:,2) = IM.get_relative_index(agents);

n = size(DATA,1);
if ~exist('prior_on_K','var')    
    prior_on_K = n-1;
end

if ~exist('gmm','var')
    gmm = gmmvar(DATA(:,1),prior_on_K);
end

centroids = cat(2,gmm.post.Norm_Mu);
Y = gmm.pjgx;


Y = Y';
% for i=1:n
%     Y(:,i) = Y(:,i) ./ sum(Y(:,i));
% end

Y_hard = get_hard_incidence_matrix(Y);

Y_hard = Y_hard';
non_active_clusters = sum(Y_hard)==0;
Y_hard(:,non_active_clusters) = [];
centroids(non_active_clusters) = [];

EVENT_TIMES = get_T_matrix(Y_hard,DATA);
Y_hard = sparse(Y_hard');

B_hard_incidence_matrix = build_B_from_Y(DATA,Y_hard);
A_hard_cooccurences = get_coocurences_in_bipartite_graph(B_hard_incidence_matrix);

%
output = struct('centroids',centroids...
    ,'A_hard_cooccurences',A_hard_cooccurences,'B_hard_incidence_matrix',B_hard_incidence_matrix,...
    'EVENT_TIMES',EVENT_TIMES,'Y_hard',Y_hard);
end
