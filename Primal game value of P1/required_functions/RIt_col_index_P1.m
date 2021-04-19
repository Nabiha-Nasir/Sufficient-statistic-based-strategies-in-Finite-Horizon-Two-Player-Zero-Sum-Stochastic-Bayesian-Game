function [col_index] = RIt_col_index_P1(t,I,A,B,k,a,n_is1)
%%% Find the column index of R_{I_t} for specific information set I
    if t==1
        col_index=(I-1)*A+a; %at t=1, I=k1. so if k1=2 then it will start after (k1-1)*A number of columns. Here a is the currect action according to the realization plan
    else
        col_index=sum(n_is1(1:t-1))*A; %how many columns are there before this stage has started
        for s=1:t-1   %if t=5 then I=k1 a1 b1 k2 a2 b2 k3 a3 b3 k4 a4 b4 k5. In this loop first ks=k1,as=a1,bs=b1 then in second look ks=k2,as=a2,bs=b2...
            ks=I(3*s-2);
            as=I(3*s-1);
            bs=I(3*s);
            col_index=col_index+A^2*(ks-1)*B*k*(A*B*k)^(t-1-s)+A*(as-1)*B*k*(A*B*k)^(t-1-s)+A*(bs-1)*k*(A*B*k)^(t-1-s);
        end
        kt=I(end);
        col_index=col_index+(kt-1)*A+a;
    end
end

