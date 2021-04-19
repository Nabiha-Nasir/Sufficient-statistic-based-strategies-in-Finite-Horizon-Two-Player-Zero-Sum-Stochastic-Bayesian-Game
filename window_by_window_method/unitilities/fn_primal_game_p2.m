function [tau,mu] = fn_primal_game_p2(T,A,B,k,l,lm,P,Q,p,q,G)
%%% to get player 2's optimal strategy and mu
    warning('off')
    %information set of player 1 and 2
    [is1,n_is1]=info_I(T,A,B,k);
    [is2,n_is2]=info_J(T,A,B,l);

    %Contraint S_{J_t}(b_t)=Q_{a_{t-1},b_{t-1}}(l_{t-1},l_t)S_{J_{t-1}}(b_t-1)
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

    %for inequality contraints
    %Creating the equation as matrix multiplication [Ain]*[variables]=[bin]. In 
    %this equation the variables are S_{J_t},(Z_{I}_{t}},Z_{\mathcal{I}_{t+1}})
    [lset,~] = Lset(T,l);

    Ain=zeros(sum(n_is1)*A,(sum(n_is2)*B+sum(n_is1)));
    bin=zeros(sum(n_is1)*A,1);

    row_index=0;
    for t=1:T
    for i=1: length (is1{t})
        for a=1:A
            row_index=row_index+1;
            for ln=1:length(lset{t})
                for b=1:B
                    [J] = construct_Is(i,ln,t,lset,is1);
                    [col_index_SJt] = SJt_col_index_P2(t,J,A,B,l,b,n_is2);
                    Ain(row_index,col_index_SJt)=lm^(t-1)*G{(is1{t}(i,(3*t-2))),J(3*t-2)}(a,b);
                end
            end
            if t<T
                for b1=1:B
                    for kplus=1:k  
                        Iplus=[is1{t}(i,:) a b1 kplus]; %problem in this line. do I need to use different row for different b1?
                        [col_index_Isplus] = I_col_index(n_is1,A,B,k,Iplus,n_is2,t+1);
                         Ain(row_index,col_index_Isplus)=P{a,b1}(is1{t}(i,end),kplus);
                    end
                end     
            end
            Ipre=is1{t}(i,:);
            [col_index_Ipre] = I_col_index(n_is1,A,B,k,Ipre,n_is2,t);
             Ain(row_index,col_index_Ipre)=-1;
        end
    end
    end

    %Construct the objective function
    f=[zeros(1,sum(n_is2)*B) p zeros(1,sum(n_is1(2:end)))]; %zeros(1,sum(n_is2)*B) is for S_{J_t} as there is no S_{J_t} in the objective function. zeros(1,sum(n_is1(2:end))) is for Z_{I_{2:T}} as in objective function there is only Z_{I_1}

    %Rearrange every coefficient to use in linprog
    f1=[zeros(1,size(Ain,2)-size(f,2)) f ]; % To ensure that matrix of objective function and Ain,Aeq are of same size we add zeroes 
    Aeq1=padarray(Aeq,[0, (size(Ain,2)-size(Aeq,2))],0,'post'); % to use linprog Ain and Aeq must have same column number as column number indicates variables. Both equation must have same variables
    lb=[zeros(sum(n_is2)*B,1);-Inf((size(Ain,1)-sum(n_is2)*B),1)];
    ub=+Inf;
    options = optimoptions('linprog','Display','none');
    [y,v2]=linprog(f1,Ain,bin,Aeq1,beq,lb,ub,options); %y=[S_{J_t} Z_{I_t}]

    %Finding the optimal strategy from the realization plan we got from previous linprog result 
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

    %to find mu
    %to find ZI1
    for k_present=1:k
        t=1;
        info=k_present;
        [col_index_ZI1] = I_col_index(n_is1,A,B,k,info,n_is2,t);
        ZI1=y(col_index_ZI1,1);
        mu(k_present)=[ZI1];
    end
    mu=-mu; %as mu=-U_{J_1}
end

