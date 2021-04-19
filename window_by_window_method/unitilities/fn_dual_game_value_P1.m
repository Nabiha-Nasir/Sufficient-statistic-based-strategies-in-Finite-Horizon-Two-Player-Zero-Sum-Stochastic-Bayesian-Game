function [Ain2,bin2,Aeq1,beq,lb,ub] = fn_dual_game_value_P1(T,A,B,k,l,lm,G,P,Q,p)

%information set of player 1 and 2
[is1,n_is1]=info_I(T,A,B,k);
[is2,n_is2]=info_J(T,A,B,l);

%Contraint R_{I_t}(a_t)=P_{a_{t-1},b_{t-1}}(k_{t-1},k_t)R_{I_{t-1}}(a_t-1)
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


%for inequality contraints
[kset,n_kset] = Kset(T,k);
[J,n_J]=info_J(T,A,B,l);

Ain=zeros(sum(n_is2)*B,(sum(n_is1)*A+sum(n_is2)+1)); %In column +1 is for U_0
bin=zeros(sum(n_is2)*B,1);

row_index=0;
for t=1:T
for j=1: length (J{t}) %for different J_t and different b_t there will be different row
    for b=1:B
        row_index=row_index+1;
        for kn=1:length(kset{t})
            for a=1:A
                [I] = construct_Is(j,kn,t,kset,J);
                [col_index] = RIt_col_index_P1(t,I,A,B,k,a,n_is1);
                Ain(row_index,col_index)=lm^(t-1)*G{I(3*t-2),(J{t}(j,(3*t-2)))}(a,b); %R_{I_s} coefficient
            end
        end
        %To find the term for U_{J_{s+1}}
        if t<T  %As for U_{J_{N+1}}=0
            for a1=1:A
                for lplus=1:l  
                    Jplus=[J{t}(j,:) a1 b lplus];
                    [col_index_Jplus] = J_col_index(n_is1,A,B,l,Jplus,n_is2,t+1);
                     Ain(row_index,col_index_Jplus)=Q{a1,b}(J{t}(j,end),lplus);%The coefficient for U_{J_{s+1}}
                end
            end      
        end
        %To find the term for U_{J_{s}}
        Jpre=J{t}(j,:);
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
    %bin1(row_index,1)=-nu(l_initial); %as in this case nu is a variable it will not be in the RHS
end


%Rearrange every coefficient to use in linprog
Ain2=[Ain;Ain1]; %Both inequality constraints need to be in one matrix thats why joining them. As linprog use <= thats why -ve sign
bin2=[bin;bin1]; %As linprog use <= thats why -ve
Aeq1=padarray(Aeq,[0, (size(Ain2,2)-size(Aeq,2))],0,'post'); % to use linprog Ain and Aeq must have same column number as column number indicates variables. Both equation must have same variables
lb=[zeros(sum(n_is1)*A,1);-Inf((size(Ain2,1)-sum(n_is1)*A),1)];
ub=+Inf;
