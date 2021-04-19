function [col_index_tau] = tau_col_index_new(A,B,l,l_present,Hb,t,n_is2)
if t==1
    col_index_tau=l_present; 
else
    col_index_tau=sum(n_is2(1:t-1));
    for s=1:t-1
        ls=Hb(3*s-2);
        as=Hb(3*s-1);
        bs=Hb(3*s);
        col_index_tau=col_index_tau+(ls-1)*A*B*l*(A*B*l)^(t-1-s)+(as-1)*B*l*(A*B*l)^(t-1-s)+(bs-1)*l*(A*B*l)^(t-1-s);
    end
    lt=Hb(end);
    col_index_tau=col_index_tau+(lt-1)+1;
end
