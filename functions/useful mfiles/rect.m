function y = rect(f)
% sets negative entries of f to 0. 

pos = (f > 0);
y = pos.*f;