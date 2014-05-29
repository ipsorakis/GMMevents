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
% total_individuals: is the total number of distinct individuals in the data stream. If not
% provided, it is automatically calculated by the second column of DATA
%
% NOTE: this code does not perform a significance test
%
% OUTPUTS:
% B: L-cell array, where L is the number of locations. Each cell B{l} is a
% N-by-K matrix, where N are the nodes and K are the groups/gathering
% events.
% T: L-cell array, where L is the number of locations. Each cell L{} is a
% K-by-2 matrix, where K is the number of gathering events. For each row k,
% the first column T(k,1) is the time when event-k started and T(k,2) the
% time where gathering event-k ended.

function [B T] = gmmevents_groups_only(DATA,total_individuals_global)

addpath(genpath('aux_scripts/'))

unique_indices = unique(DATA(:,2));

if ~exist('total_individuals_global','var')
    total_individuals_global = length(unique_indices);
end

total_individuals_current_DATA = length(unique_indices);

IM_current_DATA = index_manager(DATA(:,2));
DATA(:,2) = IM_current_DATA.get_relative_index(DATA(:,2));

total_locations = max(DATA(:,3));

B = cell(total_locations,1);
T = cell(total_locations,1);

for location_index = 1: total_locations
    
    location_indices = DATA(:,3) == location_index;
    if sum(location_indices) == 0
        continue;
    end
    
    if length(unique(DATA(location_indices,2)))<2
        continue;
    end
    
    DATA_LOC = DATA(location_indices,:);
    
    [output gmm IM_LOC] = infer_graph_from_datastream_mmVB(DATA_LOC);        
    
    B_LOC = output.B_hard_incidence_matrix;
    gathering_events_LOC = size(B_LOC,1);
    
    B_aux = zeros(gathering_events_LOC,total_individuals_current_DATA);
    B_aux(:,IM_LOC.original_indices) = B_LOC;
    
    B_aux2 = zeros(gathering_events_LOC,total_individuals_global);
    B_aux2(:,IM_current_DATA.original_indices) = B_aux;
        
    B{location_index} = B_aux2';
    
    T{location_index} = output.EVENT_TIMES;
end

end
