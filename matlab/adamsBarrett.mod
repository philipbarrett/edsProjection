%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% AUTO-GENERATED CODE FROM DYNARE.R 
% CREATED  2016-12-03-180114 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Dynare code for the Adams-Barrett model.
% Based on Alan Sutherland's code for Devreuz-Sutherland 2008

% Variable and parameter declarations

var A1 A2 NFA Y1 Y2 C1 C2 P1 P2 P11 P22 P12 P21 X11 X22 X12 X21 Q E cg cd rb1 rb2 R_1 R_2; 

varexo ep1 ep2 ep3 ep4 zeta;

parameters BT RH rho1 rho2 rho3 rho4 rho5 rho6 rho7 rho8 rho9 rho10 rho11 rho12 rho13 rho14 rho15 rho16 eta alph p1bar p2bar sigeps1 sigeps2 sigeps3 sigeps4 sigeps5 sigeps6 sigeps7 sigeps8 sigeps9 sigeps10 sigeps11 sigeps12 sigeps13 sigeps14 sigeps15 sigeps16 theta;

% Parameter values
alph = 0.86 ;
RH = 2 ;
p1bar = 1 ;
p2bar = 1 ;
BT = 0.99 ;
rho1 = 0.96 ;
rho2 = 0 ;
rho3 = 0 ;
rho4 = 0 ;
rho5 = 0 ;
rho6 = 0.96 ;
rho7 = 0 ;
rho8 = 0 ;
rho9 = 0 ;
rho10 = 0 ;
rho11 = 0.9 ;
rho12 = 0 ;
rho13 = 0 ;
rho14 = 0 ;
rho15 = 0 ;
rho16 = 0.9 ;
sigeps1 = 1.6e-05 ;
sigeps2 = 0 ;
sigeps3 = 0 ;
sigeps4 = 0 ;
sigeps5 = 0 ;
sigeps6 = 1.6e-05 ;
sigeps7 = 0 ;
sigeps8 = 0 ;
sigeps9 = 0 ;
sigeps10 = 0 ;
sigeps11 = 8e-06 ;
sigeps12 = 0 ;
sigeps13 = 0 ;
sigeps14 = 0 ;
sigeps15 = 0 ;
sigeps16 = 8e-06 ;
eta = 3 ;
theta = 0.01 ;
dr={'rb1','rb2'};                  

parameters af1;

af1 = 0;

model; 

%NFA = exp(rb2)*NFA(-1) + exp(Y1) - exp(C1) + exp(STEADY_STATE(Y1))*( BT*exp(-theta*C1)*af1*(exp(rb1) - exp(rb2)) + zeta );
NFA = exp(rb2)*NFA(-1) + exp(Y1) - exp(C1) + exp(STEADY_STATE(Y1))*( BT*af1*(exp(rb1) - exp(rb2)) + zeta );

% Main model equations

cd = C1 - C2 - Q / RH ;
    % Has conditional expectation zero (to 1st order)
cg = (1/2)*(C1 + C2);
    % Defines the excess return

A1 = rho1 * A1(-1) + rho2 * A2(-1) + rho3 * P11(-1) + rho4 * P22(-1) + ep1 ;
A2 = rho5 * A1(-1) + rho6 * A2(-1) + rho7 * P11(-1) + rho8 * P22(-1) + ep2 ;
P11 = rho9 * A1(-1) + rho10 * A2(-1) + rho11 * P11(-1) + rho12 * P22(-1) + ep3 ;
P22 = rho13 * A1(-1) + rho14 * A2(-1) + rho15 * P11(-1) + rho16 * P22(-1) + ep4 ;
 
exp(Y1) = exp(A1) * exp(P11) / exp(P1) ;
exp(Y2) = exp(A2) * exp(P22) / exp(P2) ;

exp(A1) = exp(X11) + exp(X21) ;
exp(A2) = exp(X22) + exp(X12) ;

exp(C1) = ( alph ^ ( 1 / eta ) * exp(X11) ^ ( ( eta - 1 ) / eta ) + 
             ( 1 - alph ) ^ ( 1 / eta ) * exp(X12) ^ ( ( eta - 1 ) / eta ) ) ^ ( eta / ( eta - 1 ) ) ;
exp(C2) = ( alph ^ ( 1 / eta ) * exp(X22) ^ ( ( eta - 1 ) / eta ) + 
             ( 1 - alph ) ^ ( 1 / eta ) * exp(X21) ^ ( ( eta - 1 ) / eta ) ) ^ ( eta / ( eta - 1 ) ) ;

alph * exp(C1) / exp(X11) = ( exp(P11) / exp(P1) ) ^ eta ;
alph * exp(C2) / exp(X22) = ( exp(P22) / exp(P2) ) ^ eta ;
(1-alph) * exp(C1) / exp(X12) = ( exp(P12) / exp(P1) ) ^ eta ;
(1-alph) * exp(C2) / exp(X21) = ( exp(P21) / exp(P2) ) ^ eta ;

E = P11 - P21 ;
E = P12 - P22 ;

Q = P2 - P1 + E ;

%exp(rb1) =  1/(exp(P1)*exp(Z1(-1))); 
%exp(rb2) =  exp(E) /(exp(P2)*exp(Z2(-1)));

rb1 = R_1(-1) - P1 + P1(-1) ;
rb2 = R_2(-1) - P1 + P1(-1)  + E - E(-1) ;

% BT*( exp(-RH*C2(+1)) * exp(rb2(+1)) / exp( Q(+1) ) ) = exp(-RH*C2);
% BT*( exp(-RH*C2(+1)) * exp(rb2(+1)) / exp( Q(+1) - Q ) ) = exp(-RH*C2);
BT*exp(-theta*C2)*( exp(-RH*C2(+1)) * exp(rb1(+1)) / exp( Q(+1) - Q ) ) = exp(-RH*C2);

BT*exp(-theta*C1)*( exp(-RH*C1(+1)) * exp(rb1(+1)) ) = exp(-RH*C1);
BT*exp(-theta*C1)*( exp(-RH*C1(+1)) * exp(rb2(+1)) ) = exp(-RH*C1);

end;

% Steady state values

steady_state_model;
cd = 0;
cg = 0;
NFA = 0;
Y1 = 0;
Y2 = 0;
A1 = 0;
A2 = 0;
C1 = 0;
C2 = 0;
P1 = 0; % log(p1bar);
P2 = 0; %log(p2bar);
P11 = 0 ;
P22 = 0 ;
P12 = 0 ;
P21 = 0 ;
rb1 = log(1/BT);
rb2 = log(1/BT);
R_1 = log(1/BT);
R_2 = log(1/BT);
X11 = log(alph) ;
X22 = log(alph) ;
X12 = log(1-alph) ;
X21 = log(1-alph) ;
Q = 0 ;
E = 0 ;
end; 

steady;

% Shock variances

shocks;
var ep1 = sigeps1 ;
var ep1, ep2 = sigeps2 ;
var ep1, ep3 = sigeps3 ;
var ep1, ep4 = sigeps4 ;
var ep2 = sigeps6 ;
var ep2, ep3 = sigeps7 ;
var ep2, ep4 = sigeps8 ;
var ep3 = sigeps11 ;
var ep3, ep4 = sigeps12 ;
var ep4 = sigeps16 ;
end;


% Stage 1: calculate zero-order asset holdings
% check(qz_zero_threshold=1e-12);
stoch_simul(order=1,nomoments,noprint,irf=0) ;

BB=oo_.dr.ghu;
nvar=M_.endo_nbr;
nu=M_.exo_nbr-1;
ordo=oo_.dr.order_var;

for i=1:nvar
  varo{i,1} = M_.endo_names(ordo(i),:);
  ys(i) = oo_.dr.ys(ordo(i),:);
end

SIGMA=M_.Sigma_e(1:nu,1:nu);

icd=inx('cd',varo);

[dm,na]=size(dr);

for i=1:na
ia(i)=inx(dr{i},varo);
end

D1=BB(icd,nu+1);
D2=BB(icd,1:nu);

for k=1:na-1
R1(k,:)=BB(ia(k),nu+1)-BB(ia(na),nu+1);
R2(k,:)=BB(ia(k),1:nu)-BB(ia(na),1:nu);
end

alpha_tilde=inv(R2*SIGMA*D2'*R1'-R2*SIGMA*R2'*D1)*R2*SIGMA*D2';
                
                % Set af1, af2 etc to calculated values 
                
                npr=M_.param_nbr-na+1;
                
                for i=1:na-1
                alval{i,1}=['M_.params( ' num2str(i+npr) ' )=alpha_tilde(' num2str(i) ',1);'];
                alval{i,2}=['af' num2str(i) '=alpha_tilde(' num2str(i) ',1);'];
                end   
                
                for i=1:na-1
                eval(alval{i,1});
                eval(alval{i,2});
                end  
                
                % Stage 2: calculate portfolio (first-order) and excess return (third-order) dynamics 
                % The gamma vector
                
                stoch_simul(order=2,nomoments, noprint,irf=0);
                
                BB=oo_.dr.ghu;
                EE=oo_.dr.ghxu;
                
                nx=M_.npred;
                
                D1=BB(icd,nu+1);
                D2=BB(icd,1:nu);
                
                for i=1:nx;
                for j=1:nu
                D5(i,j)=EE(icd,(i-1)*(nu+1)+j);   
                end    
                end;
                
                for k=1:na-1
                R2(k,:)=BB(ia(k),1:nu)-BB(ia(na),1:nu);
                R2k=BB(ia(k),1:nu)-BB(ia(na),1:nu);
                for i=1:nx;
                for j=1:nu
                R5k(i,j)=EE(ia(k),(i-1)*(nu+1)+j)-EE(ia(na),(i-1)*(nu+1)+j);
                end    
                end;
                CX2(k,:)=R2k*SIGMA*D5'+D2*SIGMA*R5k'; 
                end
                
                VX2=R2*SIGMA*R2'*D1;
                
                gamma=-inv(VX2)*CX2;
                
                % The delta vector
                
                icg=inx('cg',varo);
                
                G1=BB(icg,nu+1);
                G2=BB(icg,1:nu);
                
                for i=1:nx;
                for j=1:nu
                G5(i,j)=EE(icg,(i-1)*(nu+1)+j);   
                end    
                end;
                
                rs=RH*G2*SIGMA*R2';              % second-order excess returns
                
                Gx5=G5+G1*gamma'*R2;
                
                AA=oo_.dr.ghx;
                Rb3=AA(ia(na),:);
                
                for k=1:na-1
                R1k=BB(ia(k),nu+1)-BB(ia(na),nu+1);
                R2k=BB(ia(k),1:nu)-BB(ia(na),1:nu);
                for i=1:nx;
                for j=1:nu
                R5k(i,j)=EE(ia(k),(i-1)*(nu+1)+j)-EE(ia(na),(i-1)*(nu+1)+j);
                end    
                end;
                Rx5k=R5k+R1k*gamma'*R2;
                delta(k,:)=RH*R2k*SIGMA*Gx5'+RH*G2*SIGMA*Rx5k'+rs(k)*Rb3;
                end

%%% Save output for R manipulation %%%
ghx =oo_.dr.ghx;
ghu =oo_.dr.ghu;
nvar=M_.endo_nbr;
nsk=M_.nstatic;
save( 'Rvars.mat', 'ghx', 'ghu', 'nvar', 'nsk', 'delta', 'gamma', 'nu', 'varo', 'na', 'nx', 'ys', 'alpha_tilde' )
