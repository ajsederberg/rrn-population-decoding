function G = generateOrientedGraph_noSinks(num_nodes, num_edges, seed_val)
% num_nodes is the number of nodes on the graph
% num_edges is the number of edges in the graph
% num_edges must be larger than num_nodes
% seed_val initiates the random number generator.
rng(seed_val)

% each node must have at least 1 outgoing edge. No edges are bidirectional.
% first generate a list of unique num_edges pairs

G = false(num_nodes);
remaining_ii_vals = 1:num_nodes;

for i_edge = 1:num_edges
    if i_edge <= num_nodes
        jj_val = i_edge;    % enforce no-sink condition
        % don't allow self-connections
        allowed_ii_vals = setdiff(remaining_ii_vals, jj_val);
        % randomly choose ii_val from remaining set. If all of remaining
        % allowed values are > jj_val, then try again
        
        ii_val = allowed_ii_vals(randi(length(allowed_ii_vals), 1));
        % remove the chosen ii_val
        remaining_ii_vals = setdiff(remaining_ii_vals, ii_val);
        
        if jj_val < num_nodes && all(remaining_ii_vals > jj_val)
            % put ii_val back in remaining_ii_vals
            remaining_ii_vals = union(remaining_ii_vals, ii_val);
            
            % remove ii_val from allowed set
            allowed_ii_vals = setdiff(allowed_ii_vals, ii_val);
            % choose a new ii_val
            ii_val = allowed_ii_vals(randi(length(allowed_ii_vals), 1));
            % remove the new chosen ii_val
            remaining_ii_vals = setdiff(remaining_ii_vals, ii_val);
        end

    else
        jj_val = randi(num_nodes, 1);
        % allowed ii vals are nodes that are not already connected to jj_val
        % and are not jj_val
        allowed_ii_vals = setdiff(find(~G(jj_val, :)), jj_val);
        ii_val = allowed_ii_vals(randi(length(allowed_ii_vals), 1));        

    end

    
    G(ii_val, jj_val) = true;
    
end