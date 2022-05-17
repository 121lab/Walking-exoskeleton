
function [odes] = hopf_learn(t,x,mu,alpha,e, F,F_t)

% hopf频率自适应振荡器
% W = 振荡器频率
% mu = 振荡器振幅参数, m > 0
% e = 学习速度参数, e > 0
% alpha = 振荡器收敛速度
% F = 周期力 输入 (e.g., sin(t), cos(t))
% F_t = 离散时间
%
% System is 3 odes expressed in Cartesian coordinates
% X = hopf oscillator
% Y = hopf oscillator
% W = hebbian learning rule
F = interp1(F_t, F, t);
X = x(1);
Y = x(2);
W = x(3);
dXdt = alpha*(mu-(sqrt(X^2+Y^2))^2)*X-W*Y+e*F;
dYdt = alpha*(mu-(sqrt(X^2+Y^2))^2)*Y+W*X;
dWdt = -e*F*(Y/sqrt(X^2+Y^2));
odes = [dXdt; dYdt; dWdt];
end