k11_10_40 = linspace(0.1,0.3,31);
k1s1 = [linspace(0,0,7) linspace(0,0.1,4) k11_10_40(2:end-1) linspace(0.3,0.3,70) linspace(0.3,0.4,24) linspace(0.4,0.4,36)];
nor1_10_40 = linspace(3,2.2,31);
nors1 = [linspace(2,2,7) linspace(2,3,4) nor1_10_40(2:end-1) linspace(2.2,2.2,43) linspace(2.2,1.8,28) linspace(1.8,1.8,59)];
k12_50_60 = linspace(2.3,2.6,11);
k1s2 = [linspace(1.8,1.8,7) linspace(2.2,2.2,18) linspace(2.2,3,9) linspace(3,3,6) linspace(3,2.3,6) linspace(2.3,2.3,4) k12_50_60(1:end-1) linspace(2.6,2,6) linspace(2,2,4) linspace(2,1.8,11) linspace(1.8,1.8,29) linspace(1.8,1.7,11) linspace(1.7,1.7,49)];
nors2 = [linspace(0,0,7) linspace(0.3,0.3,18) linspace(0.3,0.4,9) linspace(0.4,0.4,26) linspace(0.4,0.5,11) linspace(0.5,0.5,99)];
for t = 1:1
   if t == 1
      k1s = k1s1;
      nors = nors1;
   elseif t == 1
       k1s = k1s2;
       nors = nors2;
   end
   for z = 0:168
       
       BGfilter(t,z,k1s(z+1),nors(z+1));
       
   end
end