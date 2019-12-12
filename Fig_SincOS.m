% figure;plot(sinc(0:0.01:4*pi));
%%
clear y
for OS=3:2:11
	xs=-0.5:0.01:6;
	for i=1:numel(xs)
    	% 	OS=21;
    	Locs=linspace(-1,1,OS*2+1);
    	Locs=Locs(2:2:end);
    	Phis=exp(1i*pi*Locs*xs(i));
    	y(OS,i)=sum(Phis)/OS;
	end
end
CLRs='rgbmc';
figure;
for OS=3:2:11
h((OS+1)/2)=plot(real(y(OS,:)).',CLRs((OS-1)/2),'LineWidth',2);
hold on;
end
h(1)=plot(sinc(xs),'k','LineWidth',3);
axis([numel(xs)*-0.1 numel(xs)*1.1 -1.2 1.2]);
legend(h,{'Sinc','OS=3','OS=5','OS=7','OS=9','OS=11'},'FontSize',11);
%%
figure;plot(-0.5:0.01:0.5,sinc(-0.5:0.01:0.5),'k','LineWidth',2);
axis([-1 1 0 1.5]);
%%
A=sinc(-0.5:0.001:0.5).^1;
F=ifft1cg(A,2)/sqrt(numel(A));

figure;plot(abs(F),'k-o','LineWidth',2);
axis([numel(F)*[-0.1 1.1] -1 1.3]);
axis([numel(F)*0.5*[0.97 1.03] [-0.1 1.3] ]);
