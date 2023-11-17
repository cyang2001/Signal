function FunctionE()
  [y,~] = FunctionD;
  [x,~] = FunctionA;
  [c_1, lags_1] = MyAutoCorr(x);
  [c_2, lags_2] = xcorr(x);
  [c_3, lags_3] = MyInterCorr(x,y);
  [c_4, lags_4] = xcorr(x,y);
  figure;
  subplot(4,1,1);
  stem(c_1, lags_1);
  grid on;
  subplot(4,1,2);
  stem(c_2, lags_2);
  grid on;
  subplot(4,1,3);
  stem(lags_3,c_3);
  grid on;
  subplot(4,1,4);
  stem(lags_4,c_4);
  grid on;
end
function [autoCorr,m]= MyAutoCorr(s1)
len = length(s1);
N = len;
autoCorr = zeros(1,N-1);
for m = 0:N-1
  sum=0;
  for n = 1:N - m
    sum = sum + s1(n)*s1(n+m);
  end
  autoCorr(m+1) = sum;
end
m = 0:N-1;
end
function [interCorr,m] = MyInterCorr(s1, s2)
  lenS1 = length(s1);
  lenS2 = length(s2);
  N = max(lenS1, lenS2);
  interCorr = zeros(1,2*N-1);
  for m = -N+1:N-1
    sum = 0;
    for n = 1:N
      if (n+m > 0) && (n+m <= lenS2) && (n <= lenS1)
        sum = sum + s1(n) * s2(n+m);
      end
    end
    interCorr(m+N)=sum;
  end
  m = -N+1:N-1;
end

