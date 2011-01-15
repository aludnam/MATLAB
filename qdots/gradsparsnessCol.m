function gf = gradsparsnessCol(Hkt_T)
Hkt=Hkt_T';
bar(Hkt)
sumH_t = sum(Hkt,2);
sumH_t_sq = sum(sumH_t.^2);
sumH = sum(sumH_t);
nh = length(sumH_t); 
% gf = (-1/(sqrt(nh)-1)*((sqrt(sumH_t_sq)-sumH/sqrt(sumH_t_sq))/sumH_t_sq)*sumH_t)';
gf = (1/(sqrt(nh)-1)*(sumH_t*(sumH)-sumH_t_sq)/(sumH_t_sq^1.5))';
