%% Log probability
% ここでは平均mu, 標準偏差sigmaの正規分布
% -1をかけると物理でいうポテンシャルエネルギー
function [logn] = log_normal(x, mu, sigma)
logn = -0.5.*log(2.*pi.*sigma.^2) - (x-mu).^2/(2.*sigma.^2);
end

%% Log probabilityの導関数
function [dlogn] = d_log_normal(x, mu, sigma)
dlogn = -(x-mu)./sigma.^2;
end

%% 運動エネルギー
function [Ek] = momentum(p, tau)
Ek = p.^2./(2.*tau.^2);
end

%% 運動エネルギーの導関数
function [dEk] = d_momentum(p, tau)
dEk = p./tau.^2;
end

%% ハミルトニアン
% 運動エネルギーとポテンシャルエネルギーの和
function Hamiltonian(x, p, tau, mu, sigma)
H = momentum(p, tau) + (-1.*log_normal(x, mu, sigma));
end

%% リープ・フロッグ法を事項する関数
function [x, p] = proceed_leapflog(epsilon, x, p, tau)
x = x + (-0.5.*epsilon.*(-1.*d_momentum(p, tau)));
p = p + epsilon.*d_log_normal(x, mu, sigma);
x = x + -epsilon.*(-1.*d_momentum(p, tau));
end

%% HMCを１ステップ実行するサブルーチン
function [x_accepted] = proceed_HMC_iteration(x, tau, epsilon, T)
p = randn(0, tau, [1,1]);
p_new = p;
x_new = x;
for ni = 1:T
  [x, p] = proceed_leapflog(epsilon, x_new, p_new, tau);
end
alpha = exp(Hamiltonian(x, p, tau) - Hamiltonian(x_new, p_new, tau));
u = rand(1);
if u < alpha
  x_accepted = x_new;
else
  x_accepted = x;
end

end

%% HMCを実行する関数
function [x] = proceed_HMC(tau, epsilon, T, ite, init)
x = init;
for ni = 1:ite
  x = [x, proceed_HMC_iteration(x(ni), tau, epsilon, T)];
end
end
