N = 0; %N is the number of flips
R = 0; %number of heads
H = linspace(0,1,100); %Values for hypothesis of the coin's bias

y = H.^R.*((1-H).^(N-R));
figure();
subplot(5,3,1);
plot(H,y);
xlabel('Bias weighting for heads');
ylabel('prob(H|{data},I)');
plotIterations = [1 2 3 4 8 16 32 64 128 256 512 1024 2048 4096];

for N = 1:4096 
    result = flipCoin; %1 if heads, 0 if tails
    R = R+result;
    binomial = R*log(H) + (N-R)*log(1-H);
    m = max(binomial);
    binomial = binomial - m;
    normalizingConstant = log(sum(exp(binomial)));
    y = exp(binomial - normalizingConstant);

    compare = N==plotIterations;
    indx = find(compare);
    if any(compare)
        subplot(5,3,indx+1);
        plot(H,y)
        xlabel('Bias-weighting for heads');
        ylabel('prob(H|{data},I)');
        text(0.75,0.75,num2str(plotIterations(indx)))
    end
end



function result = flipCoin
    tmp = randi(4);
    if tmp >1
        result = 0; %tails
    else
        result = 1; %heads
    end
end
