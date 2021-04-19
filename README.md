# Sufficient-statistic-based-strategies-in-Finite-Horizon-Two-Player-Zero-Sum-Stochastic-Bayesian-Game
This repository contains all the codes requires to compute sufficient statistic based strategies in a finite horizon two-player zero-sum stochastic Bayesian game:
### Problem Setup:
Jamming problems in underwater acoustic sensor networks can be modeled as Two-Player-Zero-Sum Stochastic game for infinite horizon where both players are partially informed. There are four sensors <img src="https://render.githubusercontent.com/render/math?math=(s_1,s_2,s_3,s_4)"> which are referred as Player 1 and one jammer which is referred as Player 2. Sensors' goal is to transmit data using available two channels  <img src="https://render.githubusercontent.com/render/math?math=B_1">=[10 25] kHz and <img src="https://render.githubusercontent.com/render/math?math=B_2">=[25 40] kHz. Two sensors can work at a time and each channel can be used by one sensor. Jammer's goal is to block this data transfer. The distance from sink node to sensor <img src="https://render.githubusercontent.com/render/math?math=s_1,s_2,s_3"> and <img src="https://render.githubusercontent.com/render/math?math=s_4"> are 1 km, 5 km, 1 km and 5 km, respectively. So, the states of player 1 are [1,1] (k=1), [1,5] (k=2) and [5,5] (k=3). Here, [5,1] and [1,5] are considered as the same state. The distance of the jammer from sink nodes are 0.5 km and 2 km. So, player 2 has 2 states [0.5] (l=1) and [2] (l=2). The initial state probabilities of player 1 and player 2 are p=[0.5 0.3 0.2] and q=[0.5 0.5]. The game value is the security level of the players.

The sensors coordinate with each other to use the channels so that they can maximize the data transmission. Channel <img src="https://render.githubusercontent.com/render/math?math=B_1"> is more effective for the sensors which are far away and channel <img src="https://render.githubusercontent.com/render/math?math=B_2"> is better for the sensors close by. Though the sensors and the jammer do not know each other's position, the jammer can observe whether a channel is used by a far away sensor or a close by sensor. For each time period, the jammer can only generate noises in one channel and sensors are able to detect it. If the sensor state is [1,1], the active sensors are <img src="https://render.githubusercontent.com/render/math?math=s_1"> and <img src="https://render.githubusercontent.com/render/math?math=s_3">. The feasible actions are <img src="https://render.githubusercontent.com/render/math?math=s_1"> using channel <img src="https://render.githubusercontent.com/render/math?math=B_1"> while <img src="https://render.githubusercontent.com/render/math?math=s_3"> using channel <img src="https://render.githubusercontent.com/render/math?math=B_2"> (action 1), and <img src="https://render.githubusercontent.com/render/math?math=s_1"> using channel <img src="https://render.githubusercontent.com/render/math?math=B_2"> while <img src="https://render.githubusercontent.com/render/math?math=s_3"> using channel <img src="https://render.githubusercontent.com/render/math?math=B_1"> (action 2). Similar for state [1,5] and [5,5]. The jammer's action for blocking channel <img src="https://render.githubusercontent.com/render/math?math=B_1"> is action 1 and blocking <img src="https://render.githubusercontent.com/render/math?math=B_2"> is action 2. 

As this is a stochastic game there is transition matrices for both player's state. The total channel capacity matrix (Payoff matrix), given sensors' type k and jammer's type l are as follows. 

<img src="https://user-images.githubusercontent.com/62413691/114930449-01f73200-9e03-11eb-94a9-3ca6d8a90014.png" width="400" height="150">

For example, <img src="https://render.githubusercontent.com/render/math?math=G_{k_t=2, l_t=2}(a_t=1, b_t=1)=24.89">. Transition matrix for sensor's state (P) is

<img src="https://user-images.githubusercontent.com/62413691/114935727-877de080-9e09-11eb-9943-c68c8c4e2938.png" width="300" height="150">

For example, <img src="https://render.githubusercontent.com/render/math?math=P_{a=1,b=1}(k_{t-1}=2,k_t=2)=0.4">. The transition matrix of the jammer (Q) is

<img src="https://user-images.githubusercontent.com/62413691/114936045-ec393b00-9e09-11eb-8f57-352707afa88f.png" width="200" height="100">

For example, <img src="https://render.githubusercontent.com/render/math?math=Q_{a=1,b=1}(l_{t-1}=2,l_t=2)=0.5">.

### Goal:
The goal of this project is to design a sufficient statistic based algorithm to find the optimal stratgies of the players for finite horizon game.
