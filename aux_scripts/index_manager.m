classdef index_manager < handle
    
    properties
        relative_indices;
        original_indices;
        
        original_index_lookup_table;
        relative_index_lookup_table;
        
        N;
    end
    
    methods
        %% CONSTRUCTOR
        function obj = index_manager(original_indices,relative_indices)
            
            original_indices = unique(original_indices);
            
            obj.N = length(original_indices);
            
            if ~exist('relative_indices','var')
                relative_indices = (1:obj.N)';
            end
                        
            obj.relative_indices = relative_indices;
            obj.original_indices = original_indices;
            
            obj.original_index_lookup_table = hashtable;
            obj.relative_index_lookup_table = hashtable;
            
            for i=1:obj.N
               obj.original_index_lookup_table = put(obj.original_index_lookup_table,relative_indices(i),original_indices(i));
               obj.relative_index_lookup_table = put(obj.relative_index_lookup_table,original_indices(i),relative_indices(i));
            end
        end
        %% RECOVERY FUNCTIONS
        function relative_index = get_relative_index(obj,original_index)
            n = length(original_index);
            relative_index = zeros(1,n);
            
            for i=1:n
               relative_index(i) = get(obj.relative_index_lookup_table,original_index(i)); 
            end                        
        end
        
        function original_index = get_original_index(obj,relative_index)
           n = length(relative_index);
           original_index = zeros(1,n);
           
           for i=1:n
              original_index(i) = get(obj.original_index_lookup_table,relative_index(i)); 
           end
        end        
    end
end