function [q_present] = update_belief_li_qplus(l,q, Y_star,a,b,Q)
    Y_star(isnan(Y_star))=0;
    q_temp=q.*Y_star(b,:);
    ybar=sum(q_temp);
    if ybar<10^(-7)
        q_present=zeros(1,l);
    else
        q_present=q_temp*Q{a,b}/ybar;
    end
    