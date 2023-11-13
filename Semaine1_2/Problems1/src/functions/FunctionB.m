function FunctionB()
    [x, t] = FunctionA();
    meanValue = mean(x);
    power = mean(x.^2);
    fprintf('Mean value: %f\n', meanValue);
    fprintf('Power: %f\n', power);
end
