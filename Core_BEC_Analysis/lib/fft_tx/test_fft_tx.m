% test script for fft_tx
%to do
%-  tests for unevenly sampled data

%% EVENLY SAMPLED
%test that freq amplitude values are correct in the evenly sampled case
%---start user var---
amp_strong=10;
freq_strong=100;
amp_weak=5;
freq_weak=133;
noise_amp=0;
max_time=1e2;
samp_rate=2e4;
%---end user var---
times=linspace(0,max_time,max_time*samp_rate);
testfun=@(t) amp_strong*sin(2*pi*t*freq_strong+pi)+amp_weak*sin(2*pi*t*freq_weak+pi)+noise_amp*rand(size(t));
val=testfun(times);
out=fft_tx(times,val,'padding',10,'window','gauss','win_param',{1});
figure(5)
set(gcf,'color','w')
subplot(2,1,1)
few_cycles_time=20/max([freq_strong,freq_weak]); %dont want to show all the osc only a few cycles
few_cycles_idx=ceil(few_cycles_time*samp_rate);
plot(times(1:few_cycles_idx),val(1:few_cycles_idx),'.k')
xlabel('time (s)')
ylabel('amp')
subplot(2,1,2)
plot(out(1,:),abs(out(2,:)),'k')
xlabel('Freq (Hz)')
ylabel('Amp')
xlim([0,2*max([freq_strong,freq_weak])])
%set(gca, 'YScale', 'log') %log scale
[max_amp,nearest_idx]=max(abs(out(2,:)));
max_freq=out(1,nearest_idx);
logic_str = {'FAIL', 'pass'};
amp_tol=1e-2; %fractional
freq_tol=abs(out(1,1)-out(1,2));%one freq bin %absolute
fprintf('INFO: peak freq error      : %g\n',max_freq-freq_strong)
fprintf('INFO: peak amp  frac error : %g\n',max_amp/amp_strong-1)
fprintf('TEST: peak freq within tolerance: %s\n',logic_str{1+(abs(max_freq-freq_strong)<freq_tol)})
fprintf('TEST: peak amp  within tolerance: %s\n',logic_str{1+(abs(max_amp/amp_strong-1)<amp_tol)})

%% Randomly SAMPLED
%test that freq amplitude values are correct in the randomly sampled case
% random sampling is a cool way to get closer/below the (average) niquist limit
%---start user var---
amp_strong=10;
freq_strong=100;
amp_weak=5;
freq_weak=133;
noise_amp=0;
max_time=1e1;
samp_rate=2e3;
%---end user var---
times=max_time*rand(1,max_time*samp_rate);
testfun=@(t) amp_strong*sin(2*pi*t*freq_strong+pi)+amp_weak*sin(2*pi*t*freq_weak+pi)+noise_amp*rand(size(t));
val=testfun(times);
[~,sort_idx]=sort(times);

figure(5)
set(gcf,'color','w')
subplot(2,1,1)
few_cycles_time=20/max([freq_strong,freq_weak]); %dont want to show all the osc only a few cycles
few_cycles_idx=ceil(few_cycles_time*samp_rate);
plot(times(sort_idx(1:few_cycles_idx)),val(sort_idx(1:few_cycles_idx)),'.k')
xlabel('time (s)')
ylabel('amp')

out=fft_tx(times,val,'padding',10,'window','gauss','win_param',{1});

subplot(2,1,2)
plot(out(1,:),abs(out(2,:)),'k')
xlabel('Freq (Hz)')
ylabel('Amp')
xlim([0,2*max([freq_strong,freq_weak])])
%set(gca, 'YScale', 'log') %log scale
[max_amp,nearest_idx]=max(abs(out(2,:)));
max_freq=out(1,nearest_idx);
logic_str = {'FAIL', 'pass'};
amp_tol=1e-2; %fractional
freq_tol=abs(out(1,1)-out(1,2));%one freq bin %absolute
fprintf('INFO: peak freq error      : %g\n',max_freq-freq_strong)
fprintf('INFO: peak amp  frac error : %g\n',max_amp/amp_strong-1)
fprintf('TEST: peak freq within tolerance: %s\n',logic_str{1+(abs(max_freq-freq_strong)<freq_tol)})
fprintf('TEST: peak amp  within tolerance: %s\n',logic_str{1+(abs(max_amp/amp_strong-1)<amp_tol)})



%% Very simple usage
% a very simple use case
times=linspace(0,1e3,1e6);
testfun=@(t) 100*sin(2*pi*t*100+pi)+1*sin(2*pi*t*133+pi);
val=testfun(times);
out=fft_tx(times,val,'padding',10,'window','gauss','win_param',{5});
figure(5)
plot(out(1,:),abs(out(2,:)))