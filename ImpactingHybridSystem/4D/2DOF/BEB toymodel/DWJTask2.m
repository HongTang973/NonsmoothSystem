clear all; close all; clc
% Name: Wenjun Dai
% Student number: 2037898
% Date: 20/11/2021
% Course work Task2
% The error correcting performance of a binary convolutional code

% Generator Matrices
g1 = [1,0,1,1]; % V1 = 1+x^2+x^3
g2 = [1,1,0,1]; % V2 = 1+x+x^3
g3 = [1,0,0,1]; % V3 = 1+x^3

% Generate the input numbers
N = 400;
input = randi(2,1,N)-1;

% Encode function
[Code] = encode(g1,g2,g3,input,N);

% decode function
[Path,metric] = viterbi(Code);
BER = 1-length(find(Path == input))/(N); % test the decoder

% Calculate the BER of uncoded BPSK and coded BPSK
SNRdB = -15:1:0;
SNR = 10.^(SNRdB/10);
[BER_Before] = Before(N,SNRdB);
[BER_After] = After(N,SNRdB);

%******************************plot*****************************%
figure(1)
semilogy(SNRdB,BER_Before,'r','linewidth',1.5);  % Bit Error Rate Before Coding
hold on
semilogy(SNRdB,qfunc(sqrt(2*SNR)),'g','linewidth',1.5);  % Bit Error Rate Without Coding
hold on
semilogy(SNRdB,BER_After,'b','linewidth',1.5);  % Bit Error Rate After Coding
grid on;
legend('BER Before Decoding', 'BER Without Coding', 'BER After Decoding')
xlabel('Eb/No in dB')
ylabel('BER')
title('BER in Before Decoding, After Decoding and Without Decoding in Viterbi')
axis tight


% The performance of the convolutional code in erasure channels
[BER_BEC] = BEC(N);

%% Encode function: genarate convolutional codes
function [output] = encode(g1,g2,g3,input,N)
% Calculate convolutional code
% Calculate the first value V1
% if V1 don't add the this buffer, the value of this buffer will set to 0
x0 = g1(1) * input;
x1 = g1(2) * input;
x2 = g1(3) * input;
x3 = g1(4) * input; 
% input + value of first delay buffer
x_1 = xor(x0(2:N),x1(1:N-1));
x_1 = [x0(1),x_1,x1(N)];
% using the same way like above 
x_2 = xor(x_1(3:N+1),x2(1:N-1));
x_2 = [x_1(1:2),x_2,x2(N)];
x_3 = xor(x_2(4:N+2),x3(1:N-1));
v1 = [x_2(1:3),x_3,x3(N)];

% Calculate the second value V2
% if V2 don't add the this buffer, the value of this buffer will set to 0
x0 = g2(1) * input;
x1 = g2(2) * input;
x2 = g2(3) * input;
x3 = g2(4) * input;
% input + value of first delay buffer
x_1 = xor(x0(2:N),x1(1:N-1));
x_1 = [x0(1),x_1,x1(N)];
% using the same way like above 
x_2 = xor(x_1(3:N+1),x2(1:N-1));
x_2 = [x_1(1:2),x_2,x2(N)];
x_3 = xor(x_2(4:N+2),x3(1:N-1));
v2 = [x_2(1:3),x_3,x3(N)];

% Calculate the third value V1
% if V2 don't add the this buffer, the value of this buffer will set to 0
x0 = g3(1) * input;
x1 = g3(2) * input;
x2 = g3(3) * input;
x3 = g3(4) * input;
% input + value of first delay buffer
x_1 = xor(x0(2:N),x1(1:N-1));
x_1 = [x0(1),x_1,x1(N)];
% using the same way like above 
x_2 = xor(x_1(3:N+1),x2(1:N-1));
x_2 = [x_1(1:2),x_2,x2(N)];
x_3 = xor(x_2(4:N+2),x3(1:N-1));
v3 = [x_2(1:3),x_3,x3(N)];

% organize three outputs
output = zeros(1,size(v1,2)+size(v2,2)+size(v3,2));
output(1:3:end) = v1;
output(2:3:end) = v2;
output(3:3:end) = v3;
output = output(1:length(output)-9);
%In order to flush the register, the last three bits of the input are 000. According to the code rate, the last nine bits of the output are 0, so the actual codeword discards the last nine bits.
end

%% Viterbi Decoder
function [Path,metric] = viterbi(x)
% initialize the total metrics
a_val = 0;
b_val = 0;
c_val = 0;
d_val = 0;
e_val = 0;
f_val = 0;
g_val = 0;
h_val = 0;
% set the state diagram
% a:000; b:100; c:010; d:110; e:001; f:101; g:011; h:111
% graph: [a-a:0; a-b:1; b-c:0; b-d:1; c-e:0; c-f:1; d-g:0; d-h:1; 
%         e-a:0; e-b:1; f-c:0; f-d:1; g-e:0; g-f:1; h-g:0; h-h:1]
% a-a:0 means: state a to state a should input 0
graph = [0,0,0; 1,1,1; 0,1,0; 1,0,1; 1,0,0; 0,1,1; 1,1,0; 0,0,1;
         1,1,1; 0,0,0; 1,0,1; 0,1,0; 0,1,1; 1,0,0; 0,0,1; 1,1,0];

% if there are three 000 to decode
% find the all possible ways to initialize viterbi decoder
sizex = size(x,2);
L = sizex/3;
x = [zeros(1,9),x];
% first 0 input
a_val = a_val + distance(graph(1,:),x(1:3)); % a to a
sa_1 = 0;
b_val = a_val + distance(graph(2,:),x(1:3)); % a to b
sb_1 = 1;
% second 0 input
a_val = a_val + distance(graph(1,:),x(4:6)); % a to a
sa_2 = 0;
b_val = a_val + distance(graph(2,:),x(4:6)); % a to b
sb_2 = 1;
c_val = b_val + distance(graph(3,:),x(4:6)); % b to c 
sc_2 = 0;
d_val = b_val + distance(graph(4,:),x(4:6)); % b to d
sd_2 = 1;
% third 0 input
a_val = a_val + distance(graph(1,:),x(7:9)); % a to a
sa_3 = 0;
b_val = a_val + distance(graph(2,:),x(7:9)); % a to b
sb_3 = 1;
c_val = b_val + distance(graph(3,:),x(7:9)); % b to c 
sc_3 = 0;
d_val = b_val + distance(graph(4,:),x(7:9)); % b to d
sd_3 = 1;
e_val = c_val + distance(graph(5,:),x(7:9)); % c to e
se_3 = 0;
f_val = c_val + distance(graph(6,:),x(7:9)); % c to f
sf_3 = 1;
g_val = d_val + distance(graph(7,:),x(7:9)); % d to g
sg_3 = 0;
h_val = d_val + distance(graph(8,:),x(7:9)); % d to h
sh_3 = 1;

% organize the all possible lines routes

sa = [sa_1,sa_2,sa_3,zeros(1,L)];
sb = [sa_1,sa_2,sb_3,zeros(1,L)];
sc = [sa_1,sb_2,sc_3,zeros(1,L)];
sd = [sa_1,sb_2,sd_3,zeros(1,L)];
se = [sb_1,sc_2,se_3,zeros(1,L)];
sf = [sb_1,sc_2,sf_3,zeros(1,L)];
sg = [sb_1,sd_2,sg_3,zeros(1,L)];
sh = [sb_1,sd_2,sh_3,zeros(1,L)];

% Decode cycle
for i = 4:L+3
    % a(000) comes from a(000) and e(001), compare two ways costs
    if a_val + distance(graph(1,:),x(3*i-2:3*i)) >= e_val + distance(graph(9,:),x(3*i-2:3*i))
        tempa = se;
        a_val_t = e_val + distance(graph(9,:),x(3*i-2:3*i));
        tempa(i) = 0;
    else
        tempa = sa;
        a_val_t = a_val + distance(graph(1,:),x(3*i-2:3*i));
        tempa(i) = 0;
    end
    % b(100) comes from a(000) and e(001)
    if a_val + distance(graph(2,:),x(3*i-2:3*i)) >= e_val + distance(graph(10,:),x(3*i-2:3*i))
        tempb = se;
        b_val_t = e_val + distance(graph(10,:),x(3*i-2:3*i));
        tempb(i) = 1;
    else
        tempb = sa;
        b_val_t = a_val + distance(graph(2,:),x(3*i-2:3*i));
        tempb(i) = 1;
    end
    % c(010) comes from b(100) and f(101)
    if b_val + distance(graph(3,:),x(3*i-2:3*i)) >= f_val +distance(graph(11,:),x(3*i-2:3*i))
        tempc = sf;
        c_val_t = f_val +distance(graph(11,:),x(3*i-2:3*i));
        tempc(i) = 0;
    else
        tempc = sb;
        c_val_t = b_val +distance(graph(3,:),x(3*i-2:3*i));
        tempc(i) = 0;
    end
    % d(110) comes from b(100) and f(101)
    if b_val + distance(graph(4,:),x(3*i-2:3*i)) >= f_val +distance(graph(12,:),x(3*i-2:3*i))
        tempd = sf;
        d_val_t = f_val +distance(graph(12,:),x(3*i-2:3*i));
        tempd(i) = 1;
    else
        tempd = sb;
        d_val_t = b_val +distance(graph(4,:),x(3*i-2:3*i));
        tempd(i) = 1;
    end
    % e(001) comes from c(010) and g(011)
    if c_val + distance(graph(5,:),x(3*i-2:3*i)) >= g_val +distance(graph(13,:),x(3*i-2:3*i))
        tempe = sg;
        e_val_t = g_val +distance(graph(13,:),x(3*i-2:3*i));
        tempe(i) = 0;
    else
        tempe = sc;
        e_val_t = c_val +distance(graph(5,:),x(3*i-2:3*i));
        tempe(i) = 0;
    end
    % f(101) comes from c(010) and g(011)
    if c_val + distance(graph(6,:),x(3*i-2:3*i)) >= g_val +distance(graph(14,:),x(3*i-2:3*i))
        tempf = sg;
        f_val_t = g_val +distance(graph(14,:),x(3*i-2:3*i));
        tempf(i) = 1;
    else
        tempf = sc;
        f_val_t = c_val +distance(graph(6,:),x(3*i-2:3*i));
        tempf(i) = 1;
    end
    % g(011) comes from d(110) and h(111)
    if d_val + distance(graph(7,:),x(3*i-2:3*i)) >= h_val +distance(graph(15,:),x(3*i-2:3*i))
        tempg = sh;
        g_val_t = h_val +distance(graph(15,:),x(3*i-2:3*i));
        tempg(i) = 0;
    else
        tempg = sd;
        g_val_t = d_val +distance(graph(7,:),x(3*i-2:3*i));
        tempg(i) = 0;
    end
    % h(111) comes from d(110) and h(111)
    if d_val + distance(graph(8,:),x(3*i-2:3*i)) >= h_val +distance(graph(16,:),x(3*i-2:3*i))
        temph = sh;
        h_val_t = h_val +distance(graph(16,:),x(3*i-2:3*i));
        temph(i) = 1;
    else
        temph = sd;
        h_val_t = d_val +distance(graph(8,:),x(3*i-2:3*i));
        temph(i) = 1;
    end

    % set the total metrics and path history
    a_val = a_val_t;
    b_val = b_val_t;
    c_val = c_val_t;
    d_val = d_val_t;
    e_val = e_val_t;
    f_val = f_val_t;
    g_val = g_val_t;
    h_val = h_val_t;
    sa = tempa;
    sb = tempb;
    sc = tempc;
    sd = tempd;
    se = tempe;
    sf = tempf;
    sg = tempg;
    sh = temph;
    % Compare the total metrics, find the shortest path
    % Compare a&b
    if a_val < b_val
        Path = sa;
        metric = a_val;
    else
        Path = sb;
        metric = b_val;
    end
    if c_val < metric
        Path = sc;
        metric = c_val;
    end
    if d_val < metric
        Path = sd;
        metric = d_val;
    end
    if e_val < metric
        Path = se;
        metric = e_val;
    end
    if f_val < metric
        Path = sf;
        metric = f_val;
    end
    if g_val < metric
        Path = sg;
        metric = g_val;
    end
    if h_val < metric
        Path = sh;
        metric = h_val;
    end
    % Delete the initiliza input(three 000)
    Path = Path(4:L+3);

end
end
%% Function to calculate Hamming distance between two values
function dis = distance(x,y)
dis = length(find(x ~= y));

end

%% BER of after decoding
function [BER_After] = After(N,SNRdB)

% Monte Carlo estimation
M = 100;
for i = 1:M
    % Generate the input numbers
    input(i,:) = randi(2,1,N)-1;
    % Generator Matrices
    g1 = [1,0,1,1]; % V1 = 1+x^2+x^3
    g2 = [1,1,0,1]; % V2 = 1+x+x^3
    g3 = [1,0,0,1]; % V3 = 1+x^3
    % Encode
    Code(i,:) = encode(g1,g2,g3,input(i,:),400);
    % Change 0,1 to -1,+1 to input AWGN channel
    input_AWGN(i,:) = Code(i,:)*2 -1;
end

Nr = 3*N; 
SNR = 10.^(SNRdB/10);  %Signal to Noise Ratio in Linear Scale
Eb = 1; %energy per bit
N0 = Eb ./SNR; % noise density

% input Codeword(after BPSK) to AWGN Channel
for i = 1:length(SNRdB)
    for j = 1:M
        noise(j,:) = sqrt(N0(i)/2) * randn(1,Nr);
        output(j,:) = input_AWGN(j,:) + noise(j,:);
        % Make hard decision
        for a = 1:Nr
            if output(j,a) > 0
                output(j,a) = 1;
            else
                output(j,a) = 0;
            end
        end
        % Decode
        [Path(j,:),metric(j,:)] = viterbi(output(j,:));
    end
    % Calculate BER
    BER_After(i) = 1-length(find(Path == input))/(N*M);

end
end
%% BER of Before decoding
function [BER_Before] = Before(N,SNRdB)

% Monte Carlo estimation
M = 100;
for i = 1:M
    % Generate the input numbers
    input(i,:) = randi(2,1,N)-1;
    % Generator Matrices
    g1 = [1,0,1,1]; % V1 = 1+x^2+x^3
    g2 = [1,1,0,1]; % V2 = 1+x+x^3
    g3 = [1,0,0,1]; % V3 = 1+x^3
    % Encode
    Code(i,:) = encode(g1,g2,g3,input(i,:),400);
    % Change 0,1 to -1,+1 to input AWGN channel
    input_AWGN(i,:) = Code(i,:)*2 -1;
end

Nr = 3*N; 
SNR = 10.^(SNRdB/10);  %Signal to Noise Ratio in Linear Scale
Eb = 3; %energy per bit
N0 = Eb ./SNR; % noise density

% input Codeword(after BPSK) to AWGN Channel
for i = 1:length(SNRdB)
    for j = 1:M
        noise(j,:) = sqrt(N0(i)/2) * randn(1,Nr);
        output(j,:) = input_AWGN(j,:) + noise(j,:);
        % Make hard decision
        for a = 1:Nr
            if output(j,a) > 0
                output(j,a) = 1;
            else
                output(j,a) = 0;
            end
        end
    end
    % Calculate BER
    BER_Before(i) = 1-length(find(output == Code))/(Nr*M);

end
end

%% The performance of the convolutional code in BEC
function [BER_BEC] = BEC(N)
Nr = 3*N;
% Monte Carlo estimation
M = 100;
for i = 1:M
    % Generate the input numbers
    input(i,:) = randi(2,1,N)-1;
    % Generator Matrices
    g1 = [1,0,1,1]; % V1 = 1+x^2+x^3
    g2 = [1,1,0,1]; % V2 = 1+x+x^3
    g3 = [1,0,0,1]; % V3 = 1+x^3
    % Encode
    Code(i,:) = encode(g1,g2,g3,input(i,:),400);
end

% BEC parameter(p) in [0,0.5];
Prob_error= linspace(0,1,51);
for i = 1:length(Prob_error)
    for j = 1:M
        % erasure the data randomly
        % using the value 0.1 to represent the delete bit
        Prob = rand(1,Nr); 
        Code_Del(j,:) = (Prob >= Prob_error(i)).*Code(j,:) + (Prob < Prob_error(i)).*0.5;
        % Decode the error convolutional code
        [Path(j,:),metric(j,:)] = viterbi(Code_Del(j,:));
    end
    % Calculate BER
    BER_BEC(i) = 1-length(find(Path == input))/(N*M);
end

%******************************plot*****************************%
figure(2)
plot(Prob_error,BER_BEC,'r','linewidth',1.5);
grid on;
legend('BER of convolutional code');
xlabel('Probability of Error')
ylabel('Bit Error Rate')
title('BER of convolutional code in BEC')
axis tight

end






