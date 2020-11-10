data=importdata('0.txt');
notdata=importdata('not0.txt');
h_0=data(:,2);
d_0=data(:,3);
p_0=data(:,1);
h_1=notdata(:,2);
d_1=notdata(:,3);
p_1=notdata(:,1);

hold on
plot(h_0,d_0,'.')
plot(h_1,d_1,'.')
a=[0 .1 .2 .3 .4];
%b=[.037 .037 .037 .037 .037];
c=[.1242 .1242 .1242 .1242 .1242];
%plot(b,a)
plot(a,c)
hold off
%ylabel('Mash distance')
%xlabel('Minimal hamming distance')
%axis([0 .4 0 .3])
%title('Mash distance vs minimal hamming distance')
