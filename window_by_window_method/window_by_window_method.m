clear all;
close all;
clc;
warning('off')
% In this code player 1 and 2 is playing the window-by-window method for
% finite horizon game.
%In this code by playing the game using optimal strategies for both player we are
%obtaining the game value

%% Initialization
load M.mat;%payoff matrix
G=M;
load P.mat %Player 1's transition matrix
load Q.mat %Player 2's transition matrix
N=36; %Game horizon
T=2; %window size. So in window by window it will be 1:2, 3:4
A=2; %Number of player 1's actions
B=2; %Number of player 2's actions
k=3; %Number of states of player 1
l=2; %Number of states of player 2
lm=0.3; %discounted value
p=[0.5 0.3 0.2]; %player 1's initial probability for state
q=[0.5 0.5];     %player 2's initial probability for state

% getting the first window strategy and initial vector payoff from primal game
[sigma,nu] = fn_primal_game_p1(T,A,B,k,l,lm,P,Q,p,q,G);
[tau,mu] = fn_primal_game_p2(T,A,B,k,l,lm,P,Q,p,q,G);

payoff_list=[];
payoff_detail_all=[];
my_loop=500; %number of experiemnts to get the average payoff

for loop_count=1:my_loop 
    
    payoff_detail=[];
    
    %choosing X_star and Y_star for first window
    X_star=sigma(:,1:k); %this means all rows and col 1: to k. This is the optimal strategy for player 1 at T=1
    Y_star=tau(:,1:l);%this means all rows and col 1: to L. This is the optimal strategy for player 2 at T=1
    p_present=p;
    q_present=q;
    
    stage_counter=0;
    payoff=0;
    
    %Choose state
    %Choose 1st stage's state of player 1: k_1
    [k_present] = choose_state(p,k);
    %Choose 1st stage's state of player 2: l_1
    [l_present] = choose_state(q,l);

    % Ha and Hb are the history information set player 1 and 2,
    % respectively
    Ha=[k_present];
    Hb=[l_present]; 
    action_p1=[];
    action_p2=[];
    
    %information set of player 1 and 2
    [is1,n_is1]=info_I(T,A,B,k);
    [is2,n_is2]=info_J(T,A,B,l);
     
    %% to choose actions of the players for first window
    window_stage_counter=0;
    for i=1:T
        stage_counter=stage_counter+1;
        window_stage_counter=window_stage_counter+1;
        %Choose player 1's strategy sigma
        Ha_temp=Ha;
        Ha_temp(length(Ha))=1; %As we want the starting col index of this sigma we put k_present=1
        [col_index_sigma_start] = sigma_col_index_new(A,B,k,1,Ha_temp,window_stage_counter,n_is1); %By setting k_present=1 and Ha as Ha_temp we are getting the first value of the column index for this specific sigma
        col_index_sigma=col_index_sigma_start+k_present-1; %getting col index for actual k_present
        
        %Choose player 2's strategy tau
        Hb_temp=Hb;
        Hb_temp(length(Hb))=1;%As we want the starting col index of this tau we put k_present=1
        [col_index_tau_start] = tau_col_index_new(A,B,l,1,Hb_temp,window_stage_counter,n_is2); %By setting l_present=1 and Hb as Hb_temp we are getting the first value of the column index for this specific sigma
        col_index_tau=col_index_tau_start+l_present-1; %getting col index for actual l_present
        
        %Choose player 1's action a_1 accoriding to the strategy sigma
        [a] = choose_action(sigma,A,col_index_sigma);
        action_p1=[action_p1 a];
        
        %Choose player 2's action b_1 accoriding to the strategy tau
        [b] = choose_action(tau,B,col_index_tau);
        action_p2=[action_p2 b];
        
        payoff=payoff+G{k_present,l_present}(a,b)*lm^(stage_counter-1);
        payoff_detail=[payoff_detail G{k_present,l_present}(a,b)*lm^(stage_counter-1)];

        % Find the mu value
        Y_star(isnan(Y_star))=0; % to convert the NAN values to zero
        Y_star=abs(Y_star);
        %to update mu for the next stage, call the dual game 2nd LP function with the current sufficient statistics
        [beta_vector] = fn_dual_game_2nd_P2(T,A,B,k,l,lm,G,P,Q,q_present,mu,Y_star);
        %to choose mu corresponding to the previous actions
        col_index_beta=(a-1)*B*k+(b-1)*k;
        mu=beta_vector(1,col_index_beta+1:col_index_beta+k);

        %to find q+
        [q_present] = update_belief_li_qplus(l,q_present,tau(:,col_index_tau_start:col_index_tau_start+l-1),a,b,Q);
        
        % Find the nu value
        X_star(isnan(X_star))=0; % to convert the NAN values to zero
        X_star=abs(X_star);
        %to update nu for the next stage, call the dual game 2nd LP function with the current sufficient statistics
        [alpha_vector] = fn_dual_game_2nd_P1(T,A,B,k,l,lm,G,P,Q,p_present,nu,X_star);
        %to choose nu corresponding to the previous actions
        col_index_alpha=(a-1)*B*l+(b-1)*l;
        nu=alpha_vector(1,col_index_alpha+1:col_index_alpha+l);
        
        %to find p+
        [p_present] = update_belief_li_pplus(k,p_present, sigma(:,col_index_sigma_start:col_index_sigma_start+k-1),a,b,P);

        % update state
        if stage_counter<N  %As after last stage there will be no new state
            %for k
            k_prob=P{a,b}(k_present,:);
            [k_present] = choose_state(k_prob,k);

            %for l
            l_prob=Q{a,b}(l_present,:);
            [l_present] = choose_state(l_prob,l);

        % update history
              Ha=[Ha a b k_present];
              Hb=[Hb a b l_present];
        else
            % update history
              Ha=[Ha a b];
              Hb=[Hb a b];
        end
    end
    %% for the second to last window
    while stage_counter<N
        window_stage_counter=0;
        
        % To find the size of last window
        if N-stage_counter<T
            T=N-stage_counter;
        end
        
        %to find new \tau*
        tau_plus = fn_full_dual_game_value_P2(T,A,B,k,l,lm,G,P,Q,q_present,mu);
        Y_star=tau_plus(:,1:l); %this means all rows and col 1: to L. here t=1

        %to find new \sigma*
        sigma_plus =fn_full_dual_game_value_P1(T,A,B,k,l,lm,G,P,Q,p_present,nu);
        X_star=sigma_plus(:,1:k); %this means all rows and col 1: to k. here t=1
        
        % History set for the window 
        Ha_window=[k_present]; %player 1's history set for the window
        Hb_window=[l_present]; %player 2's history set for the window
        
        for j=1:T
            window_stage_counter=window_stage_counter+1;
            stage_counter=stage_counter+1;
            
            %Choose player 2's action accoriding to the strategy updated tau
            Hb_temp=Hb_window;
            Hb_temp(length(Hb_window))=1;
            [col_index_tau_plus_start] = tau_col_index_new(A,B,l,1,Hb_temp,window_stage_counter,n_is2); %For t=1 and specific Hb this is my Y*. In this case no matter what the Hb is, it is not using it as t=1.
            col_index_tau_plus=col_index_tau_plus_start+l_present-1;
            [b] = choose_action(tau_plus,B,col_index_tau_plus);
            action_p2=[action_p2 b]; 

            %Choose player 1's action accoriding to the strategy updated sigma
            Ha_temp=Ha_window;
            Ha_temp(length(Ha_window))=1;
            [col_index_sigma_plus_start] = sigma_col_index_new(A,B,k,1,Ha_temp,window_stage_counter,n_is1); %For t=1 and specific Ha this is my X*. In this case no matter what the Ha is, it is not using it as t=1.
            col_index_sigma_plus=col_index_sigma_plus_start+k_present-1;
            [a] = choose_action(sigma_plus,A,col_index_sigma_plus);
            action_p1=[action_p1 a]; 
            
            %Compute the payoff of every stage
            payoff=payoff+G{k_present,l_present}(a,b)*lm^(stage_counter-1);
            payoff_detail=[payoff_detail G{k_present,l_present}(a,b)*lm^(stage_counter-1)];

            
            % to find updated mu
            Y_star(isnan(Y_star))=0; % to convert the NAN values to zero
            Y_star=abs(Y_star);
            %to update mu for the next stage, call the dual game 2nd LP function with the current sufficient statistics
            [beta_vector] = fn_dual_game_2nd_P2(T,A,B,k,l,lm,G,P,Q,q_present,mu,Y_star);
            %to choose mu corresponding to the previous actions
            col_index_beta=(a-1)*B*k+(b-1)*k;
            mu=beta_vector(1,col_index_beta+1:col_index_beta+k);

            %to find q+
            [q_present] = update_belief_li_qplus(l,q_present, tau_plus(:,col_index_tau_plus_start:col_index_tau_plus_start+l-1),a,b,Q);

            % to find updated nu
            X_star(isnan(X_star))=0; % to convert the NAN values to zero
            X_star=abs(X_star);
            %to update nu for the next stage, call the dual game 2nd LP function with the current sufficient statistics
            [alpha_vector] = fn_dual_game_2nd_P1(T,A,B,k,l,lm,G,P,Q,p_present,nu,X_star);    
            %to choose nu corresponding to the previous actions
            col_index_alpha=(a-1)*B*l+(b-1)*l;
            nu=alpha_vector(1,col_index_alpha+1:col_index_alpha+l);

            %to find p+
            [p_present] = update_belief_li_pplus(k,p_present,sigma_plus(:,col_index_sigma_plus_start:col_index_sigma_plus_start+k-1),a,b,P);

            % update state
            %for k
            if stage_counter<N
                k_prob=P{a,b}(k_present,:);
                [k_present] = choose_state(k_prob,k);

                %for l
                l_prob=Q{a,b}(l_present,:);
                [l_present] = choose_state(l_prob,l);

                % update history
                  Ha=[Ha a b k_present];
                  Hb=[Hb a b l_present];
                  Ha_window=[Ha_window a b k_present];
                  Hb_window=[Hb_window a b l_present];
            else
                % update history
                  Ha=[Ha a b];
                  Hb=[Hb a b];
                  Ha_window=[Ha_window a b];
                  Hb_window=[Hb_window a b];
            end  
        end
    end
    payoff_list=[payoff_list payoff];
    payoff_detail_all=[payoff_detail_all; payoff_detail];
end

payoff;
avg_payoff=sum(payoff_list)/my_loop



