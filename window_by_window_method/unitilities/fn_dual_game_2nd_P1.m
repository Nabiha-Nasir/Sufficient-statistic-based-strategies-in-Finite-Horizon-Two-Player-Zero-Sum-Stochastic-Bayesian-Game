function [alpha_vector] = fn_dual_game_2nd_P1(T,A,B,k,l,lm,G,P,Q,p,nu,X_star)

%information set of player 1 and 2
[is1,n_is1]=info_I(T,A,B,k);
[is2,n_is2]=info_J(T,A,B,l);

%As all the constraints except one are same as the previous dual value LP so we are importing previous Aeq,Ain,beq,bin from that function
[Ain_each,bin_each,Aeq_each,beq_each,lb_each,ub_each] =fn_dual_game_value_P1(T-1,A,B,k,l,lm,G,P,Q,p); %Here T-1 as the range for those equation is 1,...n-1


%As here we have (a,b)=(1,1),(1,2),(2,1),(2,2) so we need to join 4 Ain, 4 bin, 4 Aeq, 4 beq 

%for inequality constraint from dual game value function
Ain_4=zeros((size(Ain_each,1))*(A*B), (size(Ain_each,2))*(A*B));
row=0;
n=0;

for x=1:(A*B)
    Ain_row=0;
    for i=1:size(Ain_each,1)
        row=row+1;
        Ain_row=Ain_row+1;
        Ain_col=0; 
        col=size(Ain_each,2)*n;
        for j= 1:size(Ain_each,2)
            col=col+1;
            Ain_col=Ain_col+1;
            Ain_4(row,col)=Ain_each(Ain_row,Ain_col);
        end
    end
    n=n+1;
end

%for equality constraint from dual game value
Aeq_4=zeros((size(Aeq_each,1))*(A*B), (size(Aeq_each,2))*(A*B));
row=0;
n=0;

for x=1:(A*B)
    Aeq_row=0;
    for i=1:size(Aeq_each,1)
        row=row+1;
        Aeq_row=Aeq_row+1;
        Aeq_col=0;
        col=size(Aeq_each,2)*n;
        for j= 1:size(Aeq_each,2)
            col=col+1;
            Aeq_col=Aeq_col+1;
            Aeq_4(row,col)=Aeq_each(Aeq_row,Aeq_col);
        end
    end
    n=n+1;
end

beq_4=[];
%for beq
for a1=1:A
    for b1=1:B
        pplus=update_belief_li_pplus(k,p, X_star,a1,b1,P);
        pplus(isnan(pplus))=0;
        beq_each_temp=beq_each;
        beq_each_temp(1:k,1)=pplus';
        beq_4=[beq_4;beq_each_temp];
    end
end


%As in this theorm mu is a variable for the above inequality thats why adding beta
alpha_abl=zeros((size(Ain_each,1))*(A*B), A*B*l);    
row=size(Ain_each,1)-l+1;
col=1;

while row<= size(alpha_abl,1)
    for alpha_l=1:l
        alpha_abl(row,col)=1;
        row=row+1;
        col=col+1;
    end
    row=row+size(Ain_each,1)-l;
end

Ain_final=[Ain_4 alpha_abl];
Ain_all=[Ain_final zeros(size(Ain_final,1),1)]; %Here we add one extra column with zeros which represent \rho variable
bin_all=zeros(size(Ain_all,1),1);

%For varrho <= nu(l)+G_{k,l}(a,b).....
Ain2=zeros((B*l),size(Ain_all,2)); %For row l is for l no for rows for equation nu(l)_U_{J_1}>=U_0 and for column end last +1 is for U_0 variable and last +1 is for beta_ab(k)
bin2=zeros((B*l),1);

row=0;

for l_current=1:l
    x=1;
    for b=1:B
        row=row+1;
        col_varrho=size(Ain2,2); %for \varrho
        Ain2(row,col_varrho)=-1; %for \varrho
        G_constant=0;
        for a=1:A
            const=0;
            for k_current=1:k
                const=const+lm*X_star(a,k_current)*p(k_current);
                G_constant=G_constant+G{k_current,l_current}(a,b)*X_star(a,k_current)*p(k_current); %for bin2
            end
            for l_plus=1:l %for alpha
                col_alpha=size(Ain_4,2)+(a-1)*B*l+(b-1)*l+l_plus; %for alpha
                Ain2(row,col_alpha)=-const*Q{a,b}(l_current,l_plus); %for alpha
            end
            col_U0=x*size(Ain_each,2);
            Ain2(row,col_U0)=const;
            x=x+1;
        end
        bin2(row,1)=-nu(l_current)-G_constant;
    end
end

Ain_ok=-[Ain_all; Ain2];
bin_ok=-[bin_all;bin2];

%Construct the objective function
f= [zeros(1,(size(Ain_ok,2)-1)) 1];
f1=-f; %as linprog works for only minnimize obj function thats why -ve
Aeq1=padarray(Aeq_4,[0, (size(Ain_ok,2)-size(Aeq_4,2))],0,'post');
lb=[zeros(1,sum(n_is1(T-1))*A),-Inf(1,(sum(n_is2(T-1))+1)),zeros(1,sum(n_is1(T-1))*A),-Inf(1,(sum(n_is2(T-1))+1)),zeros(1,sum(n_is1(T-1))*A),-Inf(1,(sum(n_is2(T-1))+1)),zeros(1,sum(n_is1(T-1))*A),-Inf(1,(sum(n_is2(T-1))+1)), -Inf(1,(size(Ain_ok,2)-4*(size(Ain_each,2))))];
ub=+Inf*ones(size(lb));

[x,v]=linprog(f1,Ain_ok,bin_ok,Aeq1,beq_4,lb,ub);

alpha_vector=x((size(x,1)-(A*B*l)):(size(x,1)-1))';


end

