clear all;
close all;
clc;
%%% This code provides sufficient statistics based optimal strategy of
%%% player 2

%% Initialization 
load P.mat; %transition matrix of player 1
load Q.mat; %transition matrix of player 2
load M.mat; %payoff matrix
G=M;
T=2; %Number of stages in a game
A=2; %Number of player 1's actions
B=2; %Number of player 2's actions
k=3; %Number of states of player 1
l=2; %Number of states of player 2
lm=0.3; %discounted value
p=[0.5 0.3 0.2]; %player 1's initial probability for state
q=[0.5 0.5];     %player 2's initial probability for state
%to get the mu
%For example mu(1)=0.8; %initial regret over player 1's state when k=1 and mu(2)=0.2; %initial regret over player 1's state when k=2
[~,mu] = fn_primal_game_p2(T,A,B,k,l,lm,P,Q,p,q,G);

%information set of player 1 and 2
[is1,n_is1]=info_I(T,A,B,k);
[is2,n_is2]=info_J(T,A,B,l);

%% Contraint S_{J_t}(b_t)=Q_{a_{t-1},b_{t-1}}(l_{t-1},l_t)S_{J_{t-1}}(b_t-1)
%\forall t=1,...n, \forall \mathcal{J}_t
%Creating the equation as matrix multiplication [Aeq]*[variable]=[beq]. In
%this equation the variable is S_{J_t}
Aeq=zeros(sum(n_is2),sum(n_is2)*B);
beq=zeros(sum(n_is2),1);

row_index=0;
for t=1:T
    for j=1:length(is2{t})
        row_index=row_index+1;
        for b=1:B  %As there is summation b_t in S_{J_{t}}
            J=is2{t}(j,:); %for selecting each row from the information set of player 2
            [col_index_SJt] = SJt_col_index_P2(t,J,A,B,l,b,n_is2);
            Aeq(row_index,col_index_SJt)=1;
        end
        lt=J(end);
        if t==1
           beq(row_index,1)=q(lt); %As at t=1, Q_{a_0,b_0}(l_0,l_1)=1 and S_{J_0}(b_0)=Pr(l1)
        else
            Jpre=is2{t}(j,1:(length(is2{t}(j,:))-3));
            aprev=is2{t}(j,(length(is2{t}(j,:))-2));
            bprev=is2{t}(j,(length(is2{t}(j,:))-1));
            [col_index_SJtprev] = SJt_col_index_P2(t-1,Jpre,A,B,l,bprev,n_is2);
            ltpre=J(end-3);
            Aeq(row_index,col_index_SJtprev)=-Q{aprev,bprev}(ltpre,lt);
        end
    end
end
    
%% for inequality contraints
[lset,n_lset] = Lset(T,l);

Ain=zeros(sum(n_is1)*A,(sum(n_is2)*B+sum(n_is1)+1)); % In column +1 is for V_0 variable
bin=zeros(sum(n_is1)*A,1);

row_index=0;
for t=1:T
    for i=1: length (is1{t})
        for a=1:A
            row_index=row_index+1;
            for ln=1:length(lset{t})
                for b=1:B
                    [J] = construct_Is(i,ln,t,lset,is1); %As I and J are similar so we can use same function to construct both
                    [col_index_SJt] = SJt_col_index_P2(t,J,A,B,l,b,n_is2);
                    Ain(row_index,col_index_SJt)=lm^(t-1)*G{(is1{t}(i,(3*t-2))),J(3*t-2)}(a,b);
                end
            end
             %To find the term for Z_{I_{t+1}}
            if t<T %As for Z_{I_{N+1}}=0
                for b1=1:B
                    for kplus=1:k  
                        Iplus=[is1{t}(i,:) a b1 kplus]; %problem in this line. do I need to use different row for different b1?
                        [col_index_Isplus] = I_col_index(n_is1,A,B,k,Iplus,n_is2,t+1);
                         Ain(row_index,col_index_Isplus)=P{a,b1}(is1{t}(i,end),kplus); %The coefficient for Z_{I_{t+1}}
                    end
                end     
            end
            Ipre=is1{t}(i,:);
            [col_index_Ipre] = I_col_index(n_is1,A,B,k,Ipre,n_is2,t);
             Ain(row_index,col_index_Ipre)=-1;
        end
    end
end

%for mu(k)+Z_{I_1} <= V_0
Ain1=zeros(k,(sum(n_is2)*B+sum(n_is1)+1)); %For row l is for l no for rows for equation nu(l)_U_{J_1}>=U_0 and for column +1 is for U_0 variable
bin1=zeros(k,1);

row_index=0;
col_index_I1=0;
t=1; 

for k_initial=1:k
    row_index=row_index+1;
    I1=k_initial;
    [col_index_I1] = I_col_index(n_is1,A,B,k,I1,n_is2,t);
    Ain1(row_index,col_index_I1)=1;
    Ain1(row_index,size(Ain1,2))=-1; %for V_0
    bin1(row_index,1)=-mu(k_initial);
end

%% Construct the objective function

f= [zeros(1,(sum(n_is2)*B+sum(n_is1))) 1]; %zeros(1,(sum(n_is2)*B+sum(n_is1))) is for S and Z_I as they are not in the objective function. 1 is for Z_0

%Rearrange every coefficient to use in linprog
Ain2=[Ain;Ain1]; %Both inequality constraints need to be in one matrix thats why joining them
bin2=[bin;bin1]; %Both inequality constraints need to be in one matrix thats why joining them
Aeq1=padarray(Aeq,[0, (size(Ain2,2)-size(Aeq,2))],0,'post'); % to use linprog Ain and Aeq must have same column number as column number indicates variables. Both equation must have same variables
lb=[zeros(sum(n_is2)*B,1);-Inf((size(Ain2,1)-sum(n_is2)*B),1)];
ub=+Inf;
options = optimoptions('linprog','Display','none');
[y,v2]=linprog(f,Ain2,bin2,Aeq1,beq,lb,ub,options); %y=[Z_0^{a,b} Z_{I}^{a,b} S^{a,b} \beta_{a,b} \rho] and v2 should be always zero for the special mu

%% Finding the optimal strategy from the realization plan we got from previous linprog result 
%tau=SJT/(Q_{a_{t-1},b_{t-1}}(l_{t-1},l_t)*SJ_{t-1})

tau=zeros(1,sum(n_is2));
tau_col=0;

for t=1:T
    for j=1:length(is2{t})
        tau_col=tau_col+1;
        J=is2{t}(j,:); %for selecting each row from the information set of player 2
        for b=1:B
            [row_index_SJt] = SJt_col_index_P2(t,J,A,B,l,b,n_is2);
            SJt=y(row_index_SJt,1);
        lt=J(end);
        if t==1
            SJtprev=q(lt);
            qvalue=1;
        else
            Jpre=is2{t}(j,1:(length(is2{t}(j,:))-3));
            aprev=is2{t}(j,(length(is2{t}(j,:))-2));
            bprev=is2{t}(j,(length(is2{t}(j,:))-1));
            [row_index_SJtprev] = SJt_col_index_P2(t-1,Jpre,A,B,l,bprev,n_is2);
            ltpre=J(end-3);
            SJtprev=y(row_index_SJtprev,1);
            qvalue=Q{aprev,bprev}(ltpre,lt);
        end
        tau(b,tau_col)=SJt/(SJtprev*qvalue);
       end
    end
end
%%% Player 2's optimal strategy
tau %tau(1st row) for b=1 and tau(2nd row) for b=2. Each column indicates each information set at different stages. For example. 4th column is for T=2 (l_1=1, a_1=1, b_1=1, l_2=2)

