function [tau] = fn_full_dual_game_value_P2(T,A,B,k,l,lm,G,P,Q,q,mu)

%information set of player 1 and 2
[is1,n_is1]=info_I(T,A,B,k);
[is2,n_is2]=info_J(T,A,B,l);

%Contraint S_{J_t}(b_t)=Q_{a_{t-1},b_{t-1}}(l_{t-1},l_t)S_{J_{t-1}}(b_t-1)
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
           beq(row_index,1)=q(lt); %As at t=1, Q_{a_0,b_0}(l_0,l_1)=1 and S_{J_0}(b_0)=Pr(l_{t=1})=q(l)
        else
            Jpre=is2{t}(j,1:(length(is2{t}(j,:))-3));
            aprev=is2{t}(j,(length(is2{t}(j,:))-2));
            bprev=is2{t}(j,(length(is2{t}(j,:))-1));
            [col_index_SJtprev] = SJt_col_index_P2(t-1,Jpre,A,B,l,bprev,n_is2);
            ltpre=J(end-3);
            Aeq(row_index,col_index_SJtprev)=-Q{aprev,bprev}(ltpre,lt);
            %Aeq(row_index,col_index_SJtprev)=-Q{aprev,bprev}(lt,ltpre);
        end
    end
end

    
%for inequality contraints
[lset,n_lset] = Lset(T,l);
[I,n_I]=info_I(T,A,B,k);

Ain=zeros(sum(n_is1)*A,(sum(n_is2)*B+sum(n_is1)+1)); % In column +1 is for V_0 variable
bin=zeros(sum(n_is1)*A,1);

row_index=0;
for t=1:T
for i=1: length (I{t})
    for a=1:A
        row_index=row_index+1;
        for ln=1:length(lset{t})
            for b=1:B
                [J] = construct_Is(i,ln,t,lset,I); %As I and J are similar so we can use same function to construc both
                [col_index_SJt] = SJt_col_index_P2(t,J,A,B,l,b,n_is2);
                Ain(row_index,col_index_SJt)=lm^(t-1)*G{(I{t}(i,(3*t-2))),J(3*t-2)}(a,b);
            end
        end
        if t<T
            for b1=1:B
                for kplus=1:k  
                    Iplus=[I{t}(i,:) a b1 kplus]; %problem in this line. do I need to use different row for different b1?
                    [col_index_Isplus] = I_col_index(n_is1,A,B,k,Iplus,n_is2,t+1);
                     Ain(row_index,col_index_Isplus)=P{a,b1}(I{t}(i,end),kplus);
                end
            end     
        end
        Ipre=I{t}(i,:);
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



%Construct the objective function

f= [zeros(1,(sum(n_is2)*B+sum(n_is1))) 1]; 


%Rearrange every coefficient to use in linprog

Ain2=[Ain;Ain1]; %Both inequality constraints need to be in one matrix thats why joining them
bin2=[bin;bin1];
Aeq1=padarray(Aeq,[0, (size(Ain2,2)-size(Aeq,2))],0,'post'); % to use linprog Ain and Aeq must have same column number as column number indicates variables. Both equation must have same variables
%lb=[zeros(sum(n_is2)*B,1);-Inf((size(Ain2,1)-sum(n_is2)*B),1)];
lb=[zeros(1,sum(n_is2)*B),-Inf(1,(sum(n_is1)+1))];
ub=+Inf;

[y,v]=linprog(f,Ain2,bin2,Aeq1,beq,lb,ub); % here v is the game value and y=[S; Z_I; Z_0]. The number of strategy= no of values for S, Number of information sets= number of values for Z_I, number of values for Z_)=1

%Finding the optimal strategy from the realization plan we got from previous linprog result 
%tau=SJT/(Q_{a_{t-1},b_{t-1}}(l_{t-1},l_t)*SJ_{t-1})

tau=zeros(1,sum(n_is2));
tau_col=0;

row_index=0;
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


end

