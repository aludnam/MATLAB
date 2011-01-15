function gf = gradsparsnessCol2(Hkt_T, varargin)
k=varargin{1};
t=varargin{2};

Hkt=exp(reshape(Hkt_T, k, t));
%bar(Hkt)
sumH_t = sum(Hkt,2);
sumH_t_sq = sum(sumH_t.^2);
sumH = sum(sumH_t);
nh = length(sumH_t); 
% gf = (-1/(sqrt(nh)-1)*((sqrt(sumH_t_sq)-sumH/sqrt(sumH_t_sq))/sumH_t_sq)*sumH_t)';
% gf = (1/(sqrt(nh)-1)*(sumH_t*(sumH)-sumH_t_sq)/(sumH_t_sq^1.5))';
% gf_k1 = (1/(sqrt(nh)-1)*(sumH_t*(sumH)-sumH_t_sq)/(sumH_t_sq^1.5));
% %correct for linear...
gf_k1 = (1/(sqrt(nh)-1)*(sumH_t*(sumH)-sumH_t_sq)/(sumH_t_sq^1.5));
gf_kt = repmat(gf_k1,1,t).*Hkt;
gf=reshape(gf_kt, 1, k*t);

