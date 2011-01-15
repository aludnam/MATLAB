function [B,U,W,Aest,y] = ThinICAP( x, p, Delays, Splits, Choice, ...
                                  Projection, N, Stop, MaxIter, A)
%----Thin algorithm for the Independent Component Analysis of Positive Mixtures------
%
%  Extracts 'p' of the latent independent components of the observations
% by using as criterion the low-rank approximation to a set of cumulant 
% tensors. This version works for complex signals.
%
%  Two optimizations for the criterion are possible:
%
%  Thin-SVD: performs a simultaneous optimization w.r.t. of all the 
%            columns of U, and guarantees the monotonous ascent.
%
%  Thin-QR:  performs a hierarchical optimiz. w.r.t. the columns of U.
%            Guarantees the ascent of the first non-convergent mode.
%
%-----------------------------------------------------------------------
%
%  Examples:   B = ThinICA(x);  
%              B = ThinICA(x,p, Delays);   
%              B = ThinICA(x,p, Delays, Splits, Choice, Projection);
%              B = ThinICA(x,p,[10 1 1],  0   ,   1   , 0         );
%
%-----------------------------------------------------------------------
%
% [B,U,W]=ThinICA(x,p,Delays,Splits,Choice,Projection,N,Stop,MaxIter,A);
% 
% OUTPUT ARGUMENTS
%
%  B 		  Extraction matrix (B=U'*W)
%  U          Semiunitary ICA transformation.
%  W  		  Prewhitening transformation.
%  Aest       Estimated mixing system
%  y          Estimated sources.
%
% INPUT ARGUMENTS  (only 'x' is mandatory)
%
%  x          Matrix of observations, size(x)=[No.Sensors,No.Samples].
%  p          Number of independent componentes to extract.
%  Delays     = [d1,d2,d3] 3x1-vector of non-negative elements denoting
%             the no. of time tuples for each statistic in the contrast.
%               Second order statistics: (t,t),(t,t-1)...,(t,t-d1+1) 
%               Third order statistics:  
%               (t,t,t),(t,t-1,t-1),...,(t,t-d2+1,t-d2+1) 
%               Fourth order statistics: 
%               (t,t,t,t),(t,t-1,t-1,t-1),...,(t,t-d3+1,t-d3+1,t-d3+1) 
%  Splits     No. of splits of the observations 'x'. Only relevant 
%             for nonstationary signals.
%  Choice     Selects the type of optimization:
%                 1.  TSVD-ICA (simultaneous). [default]          
%                 0.  TQR-ICA  (hierarchical).
%  Projection onto the supersymmetric manifold which contains the 
%             solutions (only recommended under special conditions).
%  N          Total no. of latent independent components in 'x'.
%  Stop       Sensivity parameter to stop the convergence (around 1E-4).
%  MaxIter    Maximum no. of iterations allowed.
%  A          Mixing system (if defined activates the debug mode).
%
%-----------------------------------------------------------------------
%
% Related Bibliography:
%
% [1] S. Cruces, A. Cichocki, L. De Lathauwer, "Thin QR and SVD factori-
%     zations for simultaneous blind signal extraction", Proc. EUSIPCO, 
%     2004.
% [2] S. Cruces, A. Cichocki, "Combining blind source extraction and 
%     joint approximate diagonalization", Proc. fourth symp. on ICA/BSS,
%     pp. 463--468, Japan, 2003.
%
%-----------------------------------------------------------------------
%
% 	(c) Sergio Cruces, Andrzej Cichocki.
% 	    http://www.personal.us.es/sergio          e-mail:    sergio@us.es
% 	    Version: 2.0                          Previous update: 01/01/2004.
% 	    Version: 3.0                          Last update:     05/06/2006.
%
%-----------------------------------------------------------------------

[M,T]=size(x);
if ~exist('p','var'),           p           =1;             end;
if ~exist('Delays','var'),      Delays      =[1 1 1];       end;
if ~exist('Splits','var'),      Splits      =0;             end;
if ~exist('Choice','var'),      Choice      =1;             end;
if ~exist('Projection','var'),  Projection  =0;             end;
if ~exist('N','var'),           N           =M;             end;
if ~exist('Stop','var'),        Stop        =1e-4;          end;
if ~exist('MaxIter','var'),     MaxIter     =200;           end;
if exist('A','var'),            DebugMode   =1;             
else                            DebugMode   =0;             end;
global IcaLabMode;              % For compatibility with IcaLab.

Splits1=Splits+1;                    % no. of data blocks after split.
dT=fix(T/Splits1);                   % length of the splitted data.
Lag2=1:Delays(1);    L2=Delays(1);
Lag3=0:(Delays(2)-1);L3=Delays(2);
Lag4=0:(Delays(3)-1);L4=Delays(3); 
Statistics=(sign(Delays).*[2 3 4]); % Indicates the selected statistics.
qm=4;                               % Maximum statistic in the Contrast.
w=[0 Statistics/max(Statistics)];   % Weights of each statistic.
WhiteNoise=0;                       % White noise option is off.

% PREWHITENING
xx=x; 
x=x-mean(x')'*ones(1,T);
Rxx=x*x'/(T-1);
[Q,D]=eig(Rxx);
[v,i]=sort(abs(diag(D)));     
vx=flipud([v,i]);  
if WhiteNoise & M>N
  vn=mean(vx(N+1:M,1)); % Estimated noise power
else       
  vn=0;                 % Estimated noise power
end
W = diag(1./sqrt(vx(1:N,1)-vn))*Q(:,vx(1:N,2))';	
W_= Q(:,vx(1:N,2))*diag(sqrt(vx(1:N,1)-vn));	

z=W*x;                  % Prewhiten the data.

% SECOND ORDER STATISTICS OF THE OBSERVATIONS

if find(Statistics==2)  % Second order statistics
  for sp=1:Splits1  
    z_=z(:,(1+(sp-1)*dT):(sp*dT)); 
    for kk=1:L2
      ind=(sp-1)*L2+kk;
      k=Lag2(kk);      
      z_t1=z_(:,(1+k):dT);   % t
      z_td=z_(:,1:(dT-k));   % t-k
      aux=(z_t1*z_td.')/(dT-k-1);           
      aux=(aux+aux.')/2-(k==0)*vn; 
      vC2z(:,ind)=aux(:);    % Estimation.
    end
  end
end % if C2z

% INITIALIZATION 

it=0;
U=eye(N,p); 
for i=1:qm
  UU(:,:,i)=U;
  yy(:,:,i)=U'*z;
end;

% MAIN LOOP

NotConvergence=1;
MinIter=50;
while (it<MinIter | (it<MaxIter & NotConvergence))
  it=it+1; 
  
  Ar=W_*U;   
  Apos=max(Ar,0);
  Aneg=min(Ar,0); 
  Lneg=length(find(Aneg<0));
  % Gradient of the penalty function: f :=-trace{Aneg'*Aneg}
  dpenalty=-2*W_'*Aneg;                 
  
  
  
  U_old=UU(:,:,1+mod(it,qm));  
  
  % ESTIMATE 3rd AND 4rd ORDER STATISTICS
  
  for sp=1:Splits1
    z_=z(:,(1+(sp-1)*dT):(sp*dT));          
    y_d1=yy(:,(1+(sp-1)*dT):(sp*dT),1+mod(it-1,qm));      
    y_d2=yy(:,(1+(sp-1)*dT):(sp*dT),1+mod(it-2,qm));      
    y_d3=yy(:,(1+(sp-1)*dT):(sp*dT),1+mod(it-3,qm));      

    if find(Statistics==3) % Third order      
      for kk=1:L3
        ind=(sp-1)*L3+kk;      
        k=Lag3(kk);
        z_t1=z_(:,(1+k):dT);                
       C3zya(:,ind,:)=z_t1*(y_d1(:,1:(dT-k)).*y_d2(:,1:(dT-k))).'/(dT-k);
       C3zyb(:,ind,:)=z_t1*(y_d1(:,1:(dT-k)).*y_d3(:,1:(dT-k))).'/(dT-k);
       C3zyc(:,ind,:)=z_t1*(y_d2(:,1:(dT-k)).*y_d3(:,1:(dT-k))).'/(dT-k);
      end % kk
    end   % if
    
    if find(Statistics==4) % Fourth order 
      for kk=1:L4
        ind=(sp-1)*L4+kk;      
        k=Lag4(kk);
        z_t1=z_(:,(1+k):dT);            
        C4zy(:,ind,:)=...
           z_t1*(y_d1(:,1:(dT-k)).*y_d2(:,1:(dT-k)).*y_d3(:,1:(dT-k))...
          -diag(sum((y_d2.*y_d3).')/(dT-k))*y_d1(:,1:(dT-k))...
          -diag(sum((y_d1.*y_d3).')/(dT-k))*y_d2(:,1:(dT-k))...
          -diag(sum((y_d1.*y_d2).')/(dT-k))*y_d3(:,1:(dT-k))).'/(dT-k);              
      end % kk 
    end   % if
  end     % sp
    
  % Form the weighted statistics   
  
  for i=1:p
    
    if find(Statistics==2) % Second order  
      Ucta=kron(eye(N,N),UU(:,i,1+mod(it-1,qm)));          
      Uctb=kron(eye(N,N),UU(:,i,1+mod(it-1,qm)));             
      Uctc=kron(eye(N,N),UU(:,i,1+mod(it-1,qm)));             
      C2zya(:,:,i)=Ucta'*vC2z;
      C2zyb(:,:,i)=Uctb'*vC2z;
      C2zyc(:,:,i)=Uctc'*vC2z;
      MM(1:N,1:N,i,2)=w(2)*(C2zya(:,:,i)*C2zya(:,:,i)'+...
        C2zyb(:,:,i)*C2zyb(:,:,i)'+C2zyc(:,:,i)*C2zyc(:,:,i)');       
    end
    
    if find(Statistics==3) % Third order       
      MM(1:N,1:N,i,3)=w(3)*(C3zya(:,:,i)*C3zya(:,:,i)'+...
        C3zyb(:,:,i)*C3zyb(:,:,i)'+C3zyc(:,:,i)*C3zyc(:,:,i)');
    end
    
    if find(Statistics==4) % Fourth order 
      MM(1:N,1:N,i,4)=w(4)*C4zy(:,:,i)*C4zy(:,:,i)';
    end  
    
    Mi(:,:,i)=sum(MM(:,:,i,:),4);
    Czy(:,i)=Mi(:,:,i)*U_old(:,i);
    modes(it,i)=abs(U_old(:,i)'*Czy(:,i));  
    
  end % i=1:p
  
  % BALANCE
  eta=min(0.01+(it)/50,.5);  
       
  % ASCEND IN THE CONTRAST FUNCTION    
 
  Uprevious=U;
  if 1 %Choice==1                       % Simultaneous approach.       
      [Q_L,D,Q_R]=svd(Czy,0);        U1=Q_L*Q_R';
      [Q_L,D,Q_R]=svd(U+dpenalty,0); U2=Q_L*Q_R';           
      [Q_L,D,Q_R]=svd(U1*(1-eta)+eta*U2,0);  
      U=Q_L*Q_R';
      %TLS
      

      yest=U'*W*xx;      
      [my,t]=max(abs(yest)'); for i=1:p; signo(i)=sign(yest(i,t(i))+eps);end
      U=U*diag(signo);
      Ar=W_*U;
  else                               % Hierarchical approach
    [value,order]=sort(modes(it,:)); % Sort the modes... 
    modes(it,:)=value(p:-1:1);       % in descending order.    
    Czy(:,:)=Czy(:,order(p:-1:1));
    
    UU=UU(:,order(p:-1:1),:);        % idem for U(k),...,U(k-q+1).
    yy=yy(order(p:-1:1),:,:);        % idem for U(k),...,U(k-q+1).    
    [U,R]=qr(Czy,0);                 % Optimize using thin-QR.
  end            
    
  % PROJECTION 

  UU(:,:,1+mod(it,qm))=U;
  yy(:,:,1+mod(it,qm))=U'*z;
  if Projection,  
    for k=1:(qm-1); 
      UU(:,:,1+mod(it-k,qm))=U; 
      yy(:,:,1+mod(it-k,qm))=yy(:,:,1+mod(it,qm));
    end; 
  end
  
  % STOPING CRITERIA 
  
  NotConvergence=max(abs(U(:)-U_old(:)))>Stop;
  
  %-------- IcaLab Interrupt window ------------
  pause( 1/100 );
  go_next_step = findobj( 'Tag', 'alg_is_run' );
  if IcaLabMode & isempty(go_next_step)
    fprintf( '\nUser break.\n\n' ); break;
  end     % IcaLab Interrupt window
  
  %-------- Debug ------------------------------
  if DebugMode
      % EVALUATE THE CONTRAST FUNCTION
      Penalty(it)=norm(Aneg,'fro').^2;
      Contrast(it)=(1-eta)*sum(modes(it,:))-eta*Penalty(it)^2;;      
      
      figure(4);
      subplot(311);plot([Contrast',modes]);
      xlabel('Iterations'); ylabel('Contrast Function');      
      G=abs(Uprevious'*W*A); Gn=pinv(diag(max(G.')))*G;
      %G=abs(pinv(Apos)*A); Gn=pinv(diag(max(G.')))*G;
      Pindex(it)=(sum(sum(Gn))-p)/(p*N-p);
      subplot(312);semilogy(1:it,Pindex,'-',1:it,Penalty+2e-3,'r:');
      xlabel('Iterations'); ylabel('P_{Index}');
      subplot(313);
      error=abs([abs(U(:))-abs(U_old(:))]);
      semilogy([error,Stop*ones(p*N,1)]);axis([1 p*N Stop/1e3 max(U(:))])
      ylabel('Coefficients of U');ylabel('Convergence');
      drawnow;
      Lneg=length(find(Aneg<0));
      disp([Lneg,Penalty(it),Pindex(it)])
  end % Debug

end; % main loop.

% EVALUATE THE LAST ITERATION

it=it+1;    
for i=1:p;
  modes(it,i)=abs(U(:,i)'*Mi(:,:,i)*U(:,i));
end

[value,order]=sort(modes(it,:));   % Sort the modes... 
modes(it,:)=value(p:-1:1);         % in descending order.           
Usorted=U(:,order(p:-1:1));        % Permute the columns.

Ar=W_*Usorted; 
Apos=max(Ar,0);
Aneg=min(Ar,0);  
Lneg=length(find(Aneg<0));
if DebugMode
    % EVALUATE THE CONTRAST FUNCTION        
    Penalty(it)=norm(Aneg,'fro').^2;
    Contrast(it)=(1-eta)*sum(modes(it,:))-eta*Penalty(it);
       
    figure(4);
    subplot(311);plot([Contrast',modes]);
    xlabel('Iterations'); ylabel('Contrast Function');
    G=abs(U'*W*A); Gn=pinv(diag(max(G.')))*G;
    %Apos=max(Ar,0);
    %G=abs(pinv(Apos)*A); Gn=pinv(diag(max(G.')))*G;    
    Pindex(it)=(sum(sum(Gn))-p)/(p*N-p);    
    subplot(312);semilogy(1:it,Pindex,'-',1:it,Penalty+2e-3,'r:');
    xlabel('Iterations'); ylabel('P_{Index}');
    subplot(313);
    error=abs([abs(U(:))-abs(U_old(:))]);
    semilogy([error,Stop*ones(p*N,1)]);axis([1 p*N Stop/1e3 max(U(:))])
    ylabel('Coefficients of U');ylabel('Convergence');
    drawnow;
    disp([Lneg,Penalty(it),Pindex(it)])
end

U=Usorted;
Aest=W_*U;     % Estimated mixing system
B=Usorted'*W;  % Return the composite extraction matrix.
y=B*xx;        % Outputs

function V=unvec(A)
n=length(A);
V=reshape(A,sqrt(n),sqrt(n));


function v=vec(A)
v=A(:);
