function dc = data_circ(data,m);
% converts data to mod 2*pi format based on min distance from m

[L,D] = size(data);

if length(m) ~= D
    error('Dimension mismatch');
end;

dc = data;

for d=1:D
  e1 = abs(data(:,d) - m(d));
  e2 = abs(data(:,d)+2*pi - m(d));
  e3 = abs(data(:,d)-2*pi - m(d));
  f = find(e2 < e1);
  dc(f,d) = data(f,d) + 2*pi;
  
  f = find(e3 < e1);
  dc(f,d) = data(f,d) - 2*pi;
end;
