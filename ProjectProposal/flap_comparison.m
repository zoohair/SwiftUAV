delta_flap = [-20 -10 0 10 20]

cl0_alpha7_17p = [.316 .482 1.118 1.567 1.896];
cl0_alpha7_13p = [.409 .702 1.080 1.560 1.799];

cl0_alpha0_17p = [-.614 -.252 .295 .838 1.233];
cl0_alpha0_13p = [-.675 -.229 .325 .810 1.185];


mid = 3;
figure();
subplot(2,1,2);
plot(delta_flap,cl0_alpha7_17p - cl0_alpha7_17p(mid),'ko--','LineWidth',2); hold on;
plot(delta_flap,cl0_alpha7_13p - cl0_alpha7_13p(mid),'gx--','LineWidth',2);
title('root section at alpha = 7');
legend('t/c=17%, Re=1e6','t/c=13%, Re=3e5',4);
ylabel('delta Cl'); xlabel('delta flap');


subplot(2,1,1);
plot(delta_flap,cl0_alpha0_17p - cl0_alpha0_17p(mid),'ko--','LineWidth',2); hold on;
plot(delta_flap,cl0_alpha0_13p - cl0_alpha0_13p(mid),'gx--','LineWidth',2);
title('tip section at alpha = 0');
legend('t/c=17%, Re=1e6','t/c=13%, Re=3e5',4);
ylabel('delta Cl'); xlabel('delta flap');