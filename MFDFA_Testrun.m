fs = 10e3;
T = 10;
x = normrnd(0,1, T*fs, 1);
q = (-3:1:3)';
m = 1;

ts_min = 10/fs;
ts_max = T/6;
no_ts = 24;
ep = linspace(log2(fs*ts_min), log2(fs*ts_max), no_ts);
scale = round (2.^ep);

tic
[H1, F1] = mxMFDFA(x, scale, q, m);
toc

tic
[H2, F2] = MFDFA(x, scale, q, m);
toc