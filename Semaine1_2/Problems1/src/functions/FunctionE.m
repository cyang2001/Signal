function FunctionE()
  [y,~] = FunctionD;
  [x,~] = FunctionA;
  size(MyAutoCorr(x))
  size(xcorr(x))
  [c_3, lags_3] = MyInterCorr(x,y);
  [c_4, lags_4] = xcorr(x,y);
  figure;
  subplot(4,1,3);
  stem(lags_3,c_3);
  subplot(4,1,4);
  stem(lags_4,c_4);
end
function autoCorr = MyAutoCorr(s1)
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

