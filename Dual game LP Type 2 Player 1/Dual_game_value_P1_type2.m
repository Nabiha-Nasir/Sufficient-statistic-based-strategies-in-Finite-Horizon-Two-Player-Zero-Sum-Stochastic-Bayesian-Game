clear all;
close all;
clc;
%%% This code provides sufficient statistics based optimal strategy of
%%% player 1
%% Initialization
load M.mat;%payoff matrix
G=M;
load P.mat %Player 1's transition matrix
load Q.mat %Player 2's transition matrix
T=2; %Number of stages in a game
A=2; %Number of player 1's actions
B=2; %Number of player 2's actions
k=3; %Number of states of player 1
l=2; %Number of states of player 2
lm=0.3; %discounted value
p=[0.5 0.3 0.2]; %player 1's initial probability for state
q=[0.5 0.5];     %player 2's initial probability for state

%For example, nu(1)=0.8; %initial regret over player 2's state when l=1
%nu(2)=0.2; %initial regret over player 2's state when l=2
[~,nu] = fn_primal_game_p1(T,A,B,k,l,lm,P,Q,p,q,G);

%information set of player 1 and 2
[is1,n_is1]=info_I(T,A,B,k);
[is2,n_is2]=info_J(T,A,B,l);

%% Contraint R_{I_t}(a_t)=P_{a_{t-1},b_{t-1}}(k_{t-1},k_t)R_{I_{t-1}}(a_t-1)
Aeq=zeros(sum(n_is1),sum(n_is1)*A);
beq=zeros(sum(n_is1),1);

row_index=0;
for t=1:T
    for i=1:length(is1{t})
        row_index=row_index+1;
        for a=1:A  %As there is summation a_t in R_{I_{t}}
            I=is1{t}(i,:); %for selecting each row from the information set of player 1
            [col_index_RIt] = RIt_col_index_P1(t,I,A,B,k,a,n_is1); %finding the column index for that particular information set row
            Aeq(row_index,col_index_RIt)=1;
        end
        kt=I(end);
        if t==1
           beq(row_index,1)=p(kt); %As at t=1, P_{a_0,b_0}(k_0,k_1)=1 and R_{I_0}(a_0)=Pr(k1)
        else
            Ipre=is1{t}(i,1:(length(is1{t}(i,:))-3));
            aprev=is1{t}(i,(length(is1{t}(i,:))-2));
            bprev=is1{t}(i,(length(is1{t}(i,:))-1));
            [col_index_RItprev] = RIt_col_index_P1(t-1,Ipre,A,B,k,aprev,n_is1);
            ktpre=I(end-3);
            Aeq(row_index,col_index_RItprev)=-P{aprev,bprev}(ktpre,kt);
        end     
    end
end  


%% for inequality contraints
[kset,n_kset] = Kset(T,k);

Ain=zeros(sum(n_is2)*B,(sum(n_is1)*A+sum(n_is2)+1)); %In column +1 is for U_0
bin=zeros(sum(n_is2)*B,1);

row_index=0;
for t=1:T
    for j=1: length (is2{t}) %for different J_t and different b_t there will be different row
        for b=1:B
            row_index=row_index+1;
            for kn=1:length(kset{t})
                for a=1:A
                    [I] = construct_Is(j,kn,t,kset,is2);
                    [col_index] = RIt_col_index_P1(t,I,A,B,k,a,n_is1);
                    Ain(row_index,col_index)=lm^(t-1)*G{I(3*t-2),(is2{t}(j,(3*t-2)))}(a,b); %R_{I_t} coefficient
                end
            end
            %To find the term for U_{J_{t+1}}
            if t<T  %As for U_{J_{N+1}}=0
                for a1=1:A
                    for lplus=1:l  
                        Jplus=[is2{t}(j,:) a1 b lplus];
                        [col_index_Jplus] = J_col_index(n_is1,A,B,l,Jplus,n_is2,t+1);
                         Ain(row_index,col_index_Jplus)=Q{a1,b}(is2{t}(j,end),lplus);%The coefficient for U_{J_{t+1}}
                    end
                end      
            end
            %To find the term for U_{J_{t}}
            Jpre=is2{t}(j,:);
            [col_index_Jpre] = J_col_index(n_is1,A,B,l,Jpre,n_is2,t);
             Ain(row_index,col_index_Jpre)=-1;
        end
    end
end

%for nu(l)+U_{J_1} >= U_0
Ain1=zeros(l,(sum(n_is1)*A+sum(n_is2)+1)); %For row l is for l no for rows for equation nu(l)_U_{J_1}>=U_0 and for column +1 is for U_0 variable
bin1=zeros(l,1);

row_index=0;
col_index_J1=0;
t=1; 

for l_initial=1:l
    row_index=row_index+1;
    J1=l_initial;
    [col_index_J1] = J_col_index(n_is1,A,B,l,J1,n_is2,t);
    Ain1(row_index,col_index_J1)=1;
    Ain1(row_index,size(Ain1,2))=-1; %for U_0
    bin1(row_index,1)=-nu(l_initial);
end


%% Construct the objective function

f= [zeros(1,(sum(n_is1)*A+sum(n_is2))) 1]; %in objective function there is only \phi thats why 1 and for other variables zero

%Rearrange every coefficient to use in linprog
Ain2=-[Ain;Ain1]; %Both inequality constraints need to be in one matrix thats why joining them. As linprog use <= thats why -ve sign
bin2=-[bin;bin1]; %As linprog use <= thats why -ve
f1=-f; %as linprog works for only minnimize obj function thats why -ve
Aeq1=padarray(Aeq,[0, (size(Ain2,2)-size(Aeq,2))],0,'post'); % to use linprog Ain and Aeq must have same column number as column number indicates variables. Both equation must have same variables
lb=[zeros(sum(n_is1)*A,1);-Inf((size(Ain2,1)-sum(n_is1)*A),1)];
ub=+Inf;
options = optimoptions('linprog','Display','none');
[x,v1]=linprog(f1,Ain2,bin2,Aeq1,beq,lb,ub,options); %x=[U_0^{a,b},U_{J}^{a,b},R^{a,b},\alpha^{a,b},\phi] and v1 is the game value which should be zero for special \nu

%% Finding the optimal strategy from the realization plan we got from previous linprog result
%sigma=RIT/(P_{a_{t-1},b_{t-1}}(k_{t-1},k_t)*RI_{t-1})

sigma=zeros(1,sum(n_is1));
sigma_col=0;
for t=1:T
    for i=1:length(is1{t})
        I=is1{t}(i,:); %for selecting each row from the information set of player 1
        sigma_col=sigma_col+1;
        for a=1:A  %As there is summation a_t in R_{I_{t}}
            [row_index_RIt] = RIt_col_index_P1(t,I,A,B,k,a,n_is1); %finding the column index for that particular information set row. The number of column index in RIt is the number of row index in optimal realization plan x
            RIt=x(row_index_RIt,1);        
            kt=I(end);
            if t==1
               RItprev=p(kt);
               pvalue=1;
            else
                Ipre=is1{t}(i,1:(length(is1{t}(i,:))-3));
                aprev=is1{t}(i,(length(is1{t}(i,:))-2));
                bprev=is1{t}(i,(length(is1{t}(i,:))-1));
                [row_index_RItprev] = RIt_col_index_P1(t-1,Ipre,A,B,k,aprev,n_is1);
                ktpre=I(end-3);
                RItprev=x(row_index_RItprev,1);
                pvalue=P{aprev,bprev}(ktpre,kt);
            end
            sigma(a,sigma_col)=RIt/(RItprev*pvalue);
        end
    end
end
%%% Player 1's optimal strategy
sigma %sigma(1st row) for a=1 and sigma(2nd row) for a=2. Each column indicates each information 
%set at different stages. For example. 5th column is for T=2 (k_1=1, a_1=1, b_1=1, k_2=2). For 
%this history information set player 1 will take action a=1 with probability 0 and  action a=2 
%with probability 1

