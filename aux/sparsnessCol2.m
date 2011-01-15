function f = sparsnessCol2(Hkt_T, varargin)
k=varargin{1};
t=varargin{2};
Hkt=exp(reshape(Hkt_T, k, t));
sumH_t = sum(Hkt,2);
% figure 
% bar(sumH_t)
sumH_t_sq = sum(sumH_t.^2);
sumH = sum(sumH_t);
nh = length(sumH_t); 
f = (sqrt(nh)-sumH/(sqrt(sumH_t_sq)))/(sqrt(nh)-1);