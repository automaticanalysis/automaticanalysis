function W=sim_kt(AllD,DT,f)

MDT=trace(DT)/3;

W1111_i=Wijkl(AllD,DT,f,1,1,1,1);%1
W2222_i=Wijkl(AllD,DT,f,2,2,2,2);%2
W3333_i=Wijkl(AllD,DT,f,3,3,3,3);%3
W1112_i=Wijkl(AllD,DT,f,1,1,1,2);%4
W1113_i=Wijkl(AllD,DT,f,1,1,1,3);%5
W1222_i=Wijkl(AllD,DT,f,1,2,2,2);%6
W2223_i=Wijkl(AllD,DT,f,2,2,2,3);%7
W1333_i=Wijkl(AllD,DT,f,1,3,3,3);%8
W2333_i=Wijkl(AllD,DT,f,2,3,3,3);%9
W1122_i=Wijkl(AllD,DT,f,1,1,2,2);%10
W1133_i=Wijkl(AllD,DT,f,1,1,3,3);%11
W2233_i=Wijkl(AllD,DT,f,2,2,3,3);%12
W1123_i=Wijkl(AllD,DT,f,1,1,2,3);%13
W1223_i=Wijkl(AllD,DT,f,1,2,2,3);%14
W1233_i=Wijkl(AllD,DT,f,1,2,3,3);%15

WxD2=[  W1111_i W2222_i W3333_i W1112_i W1113_i ...
        W1222_i W2223_i W1333_i W2333_i W1122_i ...
        W1133_i W2233_i W1123_i W1223_i W1233_i];

W=WxD2/(MDT^2);