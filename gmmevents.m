% Ioannis Psorakis, Stephen J. Roberts, Iead Rezek and Ben Sheldon
% Inferring social network structure in ecological systems from spatio-temporal data streams
% University of Oxford 2011
%
% INPUTS:
% -------
% DATA: a Z-by-3 matrix, where Z the total number of observations. The
% first column Z(:,1) denotes the real-valued timestamps of each observation. The
% second column Z(:,2) denotes the indices of the individuals that were
% recorded. Such indices must be INTEGER values, as they refer to nodes in
% the adjacency matrix. Finally, the third column Z(:,3) denotes the index
% of the location where the observation took place. Again, such indices
% must be integer values.
% Example:
% % if we wish to see the contents of the 15th observation (z=15) in our
% data stream, we write:
% >> DATA(15,:)
% ans =
%    18516495           2           1
%
% which reads: "the individual i=2 appeared at location l=1 at time
% t=18516495"
%
% N: is the total number of distinct individuals in the data stream. If not
% provided, it is automatically calculated by the second column of DATA

% number_of_randomisations: the total number of data stream permutations
% during link significance testing. If the user does not wish to perform
% the test, set the value of that variable less or equal than 1, or not
% pass it at all.
%
% OUTPUTS:
% -------
% A: N-by-N adjacency matrix, where N nodes in the network
% B: L-cell array, where L is the number of locations. Each cell B{l} is a
% N-by-K matrix, where N are the nodes and K are the groups/gathering
% events.
% X: is a N-by-1 column vector, where X(n) is the number of times
% individual n appeared in the data stream
%
% Anull_mean: N-by-N adjacency matrix, were Anull_mean(i,j) is the average co-occurrence 
%(link weight) between nodes i and j under the null model
%
% Anull_std: N-by-N adjacency matrix, were Anull_std(i,j) is the stardard deviation of co-occurrence 
%(link weight) between nodes i and j under the null model
%
% Yall: L-1 cell where each element Yall{l} is a K x Z sparse matrix where each element z,k is 1 if
% observation z was assigned to cluster k and zero otherwise.
% corresponding author: Ioannis Psorakis ioannis.psorakis@eng.ox.ac.uk

function [A B X Anull_mean Anull_std Y_all] = gmmevents(DATA,total_individuals,number_of_randomisations)

addpath(genpath('aux_scripts/'))

if ~exist('total_individuals','var')
    total_individuals = length(unique(DATA(:,2)));
end

if ~exist('number_of_randomisations','var') || (exist('number_of_randomisations','var') && number_of_randomisations<=1)
    number_of_randomisations = 0;
end

total_locations = max(DATA(:,3))

A = zeros(total_individuals);

if number_of_randomisations>0
    Anull_mean = zeros(total_individuals);
    Anull_std = zeros(total_individuals);
else
    Anull_mean = nan;
    Anull_std = nan;
end

X = histc(DATA(:,2),1:total_individuals);

B = cell(total_locations,1);
Y_all = cell(total_locations,1);

for location_index = 1: total_locations
    DATA_local_worker_copy = DATA;
    location_indices = DATA_local_worker_copy(:,3) == location_index;
    if sum(location_indices) == 0
        continue;
    end
    
    if length(unique(DATA_local_worker_copy(location_indices,2)))<2
        continue;
    end
    
    DATA_LOC = DATA_local_worker_copy(location_indices,:);
    
    [output gmm IM_LOC] = infer_graph_from_datastream_mmVB(DATA_LOC);
    A_LOC = output.A_hard_cooccurences;
    B_LOC = output.B_hard_incidence_matrix;
    Y_all{location_index} = output.Y_hard
    
    %% DO SIGNIFICANCE TEST HERE
    % use either:
    if number_of_randomisations
        [A_LOC A_NULL_LOC_MEAN A_NULL_LOC_STD] = do_significance_test_of_adjancency_given_indidence_matrix(A_LOC,B_LOC,.05,number_of_randomisations);
    end
    
    %%
    A_LOC(isnan(A_LOC)) = 0;
    A(IM_LOC.original_indices,IM_LOC.original_indices) = A(IM_LOC.original_indices,IM_LOC.original_indices) + A_LOC; % or A_LOC_s
    
    B_aux = zeros(total_individuals,size(B_LOC,1));
    B_aux(IM_LOC.original_indices,:) = B_LOC';
    B{location_index} = B_aux;
    
    if number_of_randomisations
        A_NULL_LOC_MEAN(isnan(A_NULL_LOC_MEAN)) = 0;
        A_NULL_LOC_STD(isnan(A_NULL_LOC_STD)) = 0;
        Anull_mean(IM_LOC.original_indices,IM_LOC.original_indices) = Anull_mean(IM_LOC.original_indices,IM_LOC.original_indices) + A_NULL_LOC_MEAN;
        Anull_std(IM_LOC.original_indices,IM_LOC.original_indices) = Anull_std(IM_LOC.original_indices,IM_LOC.original_indices) + A_NULL_LOC_STD.^2;
    end
end

if number_of_randomisations
    Anull_std = sqrt(Anull_std);
end
end
