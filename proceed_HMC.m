%% Log probability
% �����ł͕���mu, �W���΍�sigma�̐��K���z
% -1��������ƕ����ł����|�e���V�����G�l���M�[
function [logn] = log_normal(x, mu, sigma)
logn = -0.5.*log(2.*pi.*sigma.^2) - (x-mu).^2/(2.*sigma.^2);
end

%% Log probability�̓��֐�
function [dlogn] = d_log_normal(x, mu, sigma)
dlogn = -(x-mu)./sigma.^2;
end

%% �^���G�l���M�[
function [Ek] = momentum(p, tau)
Ek = p.^2./(2.*tau.^2);
end

%% �^���G�l���M�[�̓��֐�
function [dEk] = d_momentum(p, tau)
dEk = p./tau.^2;
end

%% �n�~���g�j�A��
% �^���G�l���M�[�ƃ|�e���V�����G�l���M�[�̘a
function Hamiltonian(x, p, tau, mu, sigma)
H = momentum(p, tau) + (-1.*log_normal(x, mu, sigma));
end

%% ���[�v�E�t���b�O�@����������֐�
function [x, p] = proceed_leapflog(epsilon, x, p, tau)
x = x + (-0.5.*epsilon.*(-1.*d_momentum(p, tau)));
p = p + epsilon.*d_log_normal(x, mu, sigma);
x = x + -epsilon.*(-1.*d_momentum(p, tau));
end

%% HMC���P�X�e�b�v���s����T�u���[�`��
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

%% HMC�����s����֐�
function [x] = proceed_HMC(tau, epsilon, T, ite, init)
x = init;
for ni = 1:ite
  x = [x, proceed_HMC_iteration(x(ni), tau, epsilon, T)];
end
end
