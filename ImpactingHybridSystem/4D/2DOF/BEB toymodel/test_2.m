clear all; close all; clc

% Generate the input numbers
N = 400;% 400 bits
input = randi(2,1,N)-1;

% Generator Matrices
g1 = [1,0,1,1]; % V1 = 1+x^2+x^3
g2 = [1,1,0,1]; % V2 = 1+x+x^3
g3 = [1,0,0,1]; % V3 = 1+x^3



% Calculate convolutional code
% Calculate the first value V1
x0 = g1(1) * input;
x1 = g1(2) * input;
x2 = g1(3) * input;
x3 = g1(4) * input; 
% input + value of first delay buffer
x_1 = xor(x0(2:N),x1(1:N-1));
x_1 = [x0(1),x_1,x1(N)];
%input + value of first delay buffer+ value of second delay buffer
x_2 = xor(x_1(3:N+1),x2(1:N-1));
x_2 = [x_1(1:2),x_2,x2(N)];
%input + value of all delay buffer
x_3 = xor(x_2(4:N+2),x3(1:N-1));
v1 = [x2(1:3),x_3,x3(N)];

% Calculate the second value V2
x0 = g2(1) * input;
x1 = g2(2) * input;
x2 = g2(3) * input;
x3 = g2(4) * input;

x_1 = xor(x0(2:N),x1(1:N-1));
x_1 = [x0(1),x_1,x1(N)];
x_2 = xor(x_1(3:N+1),x2(1:N-1));
x_2 = [x_1(1:2),x_2,x2(N)];
x_3 = xor(x_2(4:N+2),x3(1:N-1));
v2 = [x_2(1:3),x_3,x3(N)];

% Calculate the third value V1
x0 = g3(1) * input;
x1 = g3(2) * input;
x2 = g3(3) * input;
x3 = g3(4) * input;

x_1 = xor(x0(2:N),x1(1:N-1));
x_1 = [x0(1),x_1,x1(N)];
x_2 = xor(x_1(3:N+1),x2(1:N-1));
x_2 = [x_1(1:2),x_2,x2(N)];
x_3 = xor(x_2(4:N+2),x3(1:N-1));
v3 = [x_2(1:3),x_3,x3(N)];

% Codeword= v1v2v3 matrix
Codeword = zeros(1,size(v1,2)+size(v2,2)+size(v3,2));
Codeword (1:3:end) = v1;
Codeword (2:3:end) = v2;
Codeword (3:3:end) = v3;
Codeword  = Codeword (1:length(Codeword ) - 9);
%%In order to flush the register, the last three bits of the input are 000. According to the code rate, the last nine bits of the output are 0, so the actual codeword discards the last nine bits.

%Decode function
[m,t] = viterbi(Codeword);
BER = 1-length(find(m == input))/(N);
sizex = size(Codeword,2);
L = sizex/3;
%% BER Before decording
M = 100; % Monte Carlo estimation of BER
for i = 1:M
input(i,:) = randi(2,1,N)-1;

 n=3;k=1;K=3;
gen=[13,15,11];
trellis= poly2trellis(L,gen);
Codeword(i,:)= convenc(input(i,:),trellis);
Codeword_BPSK(i,:) =Codeword(i,:)*2 -1; % BPSK modulation
end

Nr = 3*N; % code rate=1/3
SNRdB = -15:1:0;       %Signal to Noise Ratio in dB
SNR = 10.^(SNRdB/10);  %Signal to Noise Ratio in Linear Scale
Eb = 3;
N0 = Eb ./SNR;

for k = 1:length(SNRdB) % Monte Carlo estimation of BER
    for j = 1:M
    noise(j,:) = sqrt(N0(k)/2) * randn(1,Nr);
    output(j,:) = Codeword_BPSK(j,:) + noise(j,:);
    
    for h =  1: Nr
        if output(j,h) >= 0
            output(j,h) = 1;
        else
            output(j,h) = 0;
        end
    end

    end
    % Calculate BER
    BER_Before(k) = 1-length(find(m == input))/(N*M);
end

%% BER after decoding
M = 100; % Monte Carlo estimation of BER
for i = 1:M
input(i,:) = randi(2,1,N)-1;

 n=3;k=1;K=3;
gen=[13,15,11];
trellis= poly2trellis(L,gen);
Codeword(i,:)= convenc(input(i,:),trellis);
Codeword_BPSK(i,:) =Codeword(i,:)*2 -1; % BPSK modulation
end

Nr = 3*N; % code rate=1/3
SNRdB = -15:1:0;       %Signal to Noise Ratio in dB
SNR = 10.^(SNRdB/10);  %Signal to Noise Ratio in Linear Scale
Eb = 1;
N0 = Eb ./SNR;

for k = 1:length(SNRdB) % Monte Carlo estimation of BER
    for j = 1:M
    noise(j,:) = sqrt(N0(k)/2) * randn(1,Nr);
    output(j,:) = Codeword_BPSK(j,:) + noise(j,:);
    
    for h =  1: Nr
        if output(j,h) >= 0
            output(j,h) = 1;
        else
            output(j,h) = 0;
        end
    end
        [m(j,:),t(j,:)] = viterbi(output(j,:))
    end
    % Calculate BER
    BER_Before(k) = 1-length(find(m == Codeword))/(Nr*M);
end


%******************************plot*****************************%
semilogy(SNRdB,BER_Before,'g','linewidth',1.5);  % Bit Error Rate Before Coding
hold on
semilogy(SNRdB,qfunc(sqrt(2*SNR)),'b','linewidth',1.5);  % Bit Error Rate Without Coding
hold on
semilogy(SNRdB,BER_After,'r','linewidth',1.5);  % Bit Error Rate After Coding
grid on;
legend('BER Before Decoding', 'BER Without Coding', 'BER After Decoding')
xlabel('SNR in dB')
ylabel('BER')
title('BER in Before Decoding, After Decoding and Without Decoding in Viterbi')
axis tight


% Script function
%% Viterbi hard decision
% state a:000; b:100; c:010; d:110; e:001; f:101; g:011; h:111
function [m,t] = viterbi(x)
sizex = size(x,2);
L = sizex/3;
x = [zeros(1,9),x];
% to record the value
val_a = 0;
val_b = 0;
val_c = 0;
val_d = 0;
val_e = 0;
val_f = 0;
val_g = 0;
val_h = 0;
% graph   aa0    ab1    bc0   bd1    ce0     cf1   dg0   dh1      
graph = [0,0,0; 1,1,1; 0,1,0; 1,0,1; 1,0,0; 0,1,1; 1,1,0; 0,0,1;
         1,1,1; 0,0,0; 1,0,1; 0,1,0; 0,1,1; 1,0,0; 0,0,1; 1,1,0];
%         ea0    eb1    fc0    fd1    ge0    gf1    hhg0   hh1

% to record route 
ma = zeros(1,L+2);
mb = zeros(1,L+2);
mc = zeros(1,L+2);
md = zeros(1,L+2);
me = zeros(1,L+2);
mf = zeros(1,L+2);
mg = zeros(1,L+2);
mh = zeros(1,L+2);



%First Input 0
val_a = val_a + dis(graph(1,:),x(1:3));% a to a
ma(1)=0;% input=0
val_b = val_b + dis(graph(2,:),x(1:3));% a to b
mb(1)=1;% input=1

%Second Input 0
val_a = val_a + dis(graph(1,:),x(4:6)); % a to a
ma(2)=0;% input=0
val_b = val_b + dis(graph(2,:),x(4:6)); % a to b
mb(2)=1;% input=1
val_c = val_b + dis(graph(3,:),x(4:6)); % b to c
mc(2)=0;% input=0
val_d = val_b + dis(graph(4,:),x(4:6)); % b to d
md(2)=1;% input=1

%Third Input 0
val_a = val_a + dis(graph(1,:),x(7:9)); % a to a
ma(3) = 0;
val_b = val_a + dis(graph(2,:),x(7:9)); % a to b
mb(3) = 1;
val_c = val_b + dis(graph(3,:),x(7:9)); % b to c 
mc(3) = 0;
val_d = val_b + dis(graph(4,:),x(7:9)); % b to d
md(3) = 1;
val_e = val_c + dis(graph(5,:),x(7:9)); % c to e
me(3) = 0;
val_f = val_c + dis(graph(6,:),x(7:9)); % c to f
mf(3) = 1;
val_g = val_d + dis(graph(7,:),x(7:9)); % d to g
mg(3) = 0;
val_h = val_d + dis(graph(8,:),x(7:9)); % d to h
mh(3) = 1;

for i = 4:L+3
%     val_a_t =val_a; val_b_t =val_b;val_c_t =val_c;val_d_t =val_d;
%     val_e_t =val_e;val_f_t =val_f;val_g_t =val_g;val_h_t =val_h;
%     tempa = ma;tempb = mb;tempc = mc;tempd = md;
%     tempe = me;tempf = mf;tempg = mg;temph = mh;
% for val_a
        if val_a + dis(graph(1,:), x(3*i-2:3*i)) >= val_e + dis(graph(9,:),x(3*i-2:3*i))
            tempa = me; 
            val_a_t = val_e + dis(graph(9,:),x(3*i-2:3*i));
            tempa(i)=0;
        else
            val_a_t = val_a + dis(graph(1,:),x(3*i-2:3*i));
            tempa = ma;
            tempa(i)=0;
        end
        
%for val_b
         if val_a + dis(graph(2,:),x(3*i-2:3*i)) >= val_e + dis(graph(10,:),x(3*i-2:3*i))
            tempb = me; 
            val_b_t = val_e + dis(graph(10,:),x(3*i-2:3*i));
            tempb(i)=1;
        else
            val_b_t = val_a + dis(graph(2,:),x(3*i-2:3*i));
            tempb = ma;
            tempb(i)=1;         
         end
        
%for val_c 
            if val_b + dis(graph(3,:),x(3*i-2:3*i)) >= val_f + dis(graph(11,:),x(3*i-2:3*i))
            tempc = mf; 
            val_c_t = val_f + dis(graph(11,:),x(3*i-2:3*i));
            tempc(i)=0;
            else
            val_c_t = val_b + dis(graph(3,:),x(3*i-2:3*i));
            tempc = mb;
            tempc(i)=0;
            end
            
%for val_d
            if val_b + dis(graph(4,:),x(3*i-2:3*i)) >= val_f +dis(graph(12,:),x(3*i-2:3*i))
            tempd = mf;
            val_d_t = val_f +dis(graph(12,:),x(3*i-2:3*i));
            tempd(i) = 1;
            else
            val_d_t = val_b +dis(graph(4,:),x(3*i-2:3*i));
            tempd = mb;
            tempd(i) = 1;
            end
            
%for val_e 
            if val_c + dis(graph(5,:),x(3*i-2:3*i)) >= val_g +dis(graph(13,:),x(3*i-2:3*i))
            tempe = mg;
            val_e_t = val_g +dis(graph(13,:),x(3*i-2:3*i));
            tempe(i) = 0; 
            else
            val_e_t = val_c +dis(graph(5,:),x(3*i-2:3*i));
            tempe = mc;
            tempe(i) = 0;
            end
%for val_f 
            if val_c + dis(graph(6,:),x(3*i-2:3*i)) >= val_g +dis(graph(14,:),x(3*i-2:3*i))
            tempf = mg;
            val_f_t = val_g +dis(graph(14,:),x(3*i-2:3*i));
            tempf(i) = 1;
            else
            val_f_t = val_c+dis(graph(6,:),x(3*i-2:3*i));
            tempf = mc;
            tempf(i) = 1;
            end
%for val_g 
            if val_d + dis(graph(7,:),x(3*i-2:3*i)) >= val_h +dis(graph(15,:),x(3*i-2:3*i))
            tempg = mh;
            val_g_t = val_h +dis(graph(15,:),x(3*i-2:3*i));
            tempg(i) = 0;
            else
            val_g_t = val_d +dis(graph(7,:),x(3*i-2:3*i));
            tempg = md;
            tempg(i) = 0;
            end
            
%for val_h 
            if val_d + dis(graph(8,:),x(3*i-2:3*i)) >= val_h+dis(graph(16,:),x(3*i-2:3*i))
            temph = mh;
            val_h_t = val_h +dis(graph(16,:),x(3*i-2:3*i));
            temph(i) = 1;
            else
            val_h_t = val_d +dis(graph(8,:),x(3*i-2:3*i));
            temph = md;
            temph(i) = 1;
            end
    val_a =val_a_t;        
    val_b =val_b_t;
    val_c =val_c_t;
    val_d =val_d_t;
    ma = tempa;
    mb = tempb;
    mc = tempc;
    md = tempd;        
end

% Best path
if val_a <= val_b
    m = ma;
    t = val_a;
else
    m = mb;
   t = val_a;
end

if val_c <= t
    m = mc;
    t =val_c;
end

if val_d <= t
    m = md;
    t = val_d;
end
% Ignore the 000 used to flush registers
m = m(4:L+3);

% inner function
% function to compute distances between x,y
function d = dis(x, y)   
  d = sum(xor(x,y));
end
%
end
