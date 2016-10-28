Fe2=[0.01 0.05 0.1 0.2 0.4 0.6 0.8 0 1 2];
for i=1:1
    Se2=Fe2(i);
    for it=1:1
        %Matrix Generation
        z1=randn(10,1);
        z2=randn(10,1);
        z3=randn(10,1);
        z4=randn(10,1);
       y11=randn(10,1)*sqrt(Se2);
       y12=randn(10,1)*sqrt(Se2);
       y13=randn(10,1)*sqrt(Se2);
       y14=randn(10,1)*sqrt(Se2);
       y15=randn(10,1)*sqrt(Se2);
       y21=randn(10,1)*sqrt(Se2);
       y22=randn(10,1)*sqrt(Se2);
       y23=randn(10,1)*sqrt(Se2);
       y24=randn(10,1)*sqrt(Se2);
       y25=randn(10,1)*sqrt(Se2);
       y31=randn(10,1)*sqrt(Se2);
       y32=randn(10,1)*sqrt(Se2);
       y33=randn(10,1)*sqrt(Se2);
       y34=randn(10,1)*sqrt(Se2);
       y35=randn(10,1)*sqrt(Se2);
       y41=randn(10,1)*sqrt(Se2);
       y42=randn(10,1)*sqrt(Se2);
       y43=randn(10,1)*sqrt(Se2);
       y44=randn(10,1)*sqrt(Se2);
       y45=randn(10,1)*sqrt(Se2);
       x11=z1+y11;
       x12=z1+y12;
       x13=z1+y13;
       x14=z1+y14;
       x15=z1+y15;
       
       x21=z2+y21;
       x22=z2+y22;
       x23=z2+y23;
       x24=z2+y24;
       x25=z2+y25;
       
       x31=z2+y31;
       x32=z2+y32;
       x33=z2+y33;
       x34=z2+y34;
       x35=z2+y35;
       
        x41=z2+y41;
       x42=z2+y42;
       x43=z2+y43;
       x44=z2+y44;
       x45=z2+y45;
       
      pi= [x11 x12 x13 x14 x15 x21 x22 x23 x24 x25 x31 x32 x33 x34 x35 x41 x42 x43 x44 x45]

      
      
    end
end
       
        