function [l] = distance(x,n)
%This function takes the solution returned by the DILP.m file and converts
%the entries in the vector from pauplin coordinates to real distances
%using the formula: l_ij = n-2-log2(x_ij) for all entries
num_entries = length(x);
l = zeros(num_entries,1);
    for i = 1:num_entries
        l(i) = n-2-log2(x(i));
    end
end
