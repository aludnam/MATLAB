function vec_colum=makeveccolumn(vec)
% vec_colum=makeveccolumn(vec)
% Makes sure that the output vector vec_col is a column
sv=size(vec);
vec_colum=vec;
if sv(1)<sv(2)
    vec_colum=vec';
end