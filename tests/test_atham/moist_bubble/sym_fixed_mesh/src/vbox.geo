Point(1)={500,1500,0};
Point(2)={500,0,0};

N1=60;
N3=10;
N2=20;


num[]=Extrude {0,-1240,0} {Point{1};Layers{N1};};
num2[]=Extrude {0,260,0} {Point{2};Layers{N3};};
s11[]=Extrude {500,0,0} {Line{num[1]};Layers{N2};};
s12[]=Extrude {500,0,0} {Line{num2[1]};Layers{N2};};
s21[]=Extrude {-500,0,0} {Line{num[1]};Layers{N2};};
s22[]=Extrude {-500,0,0} {Line{num2[1]};Layers{N2};};

Physical Line(11)={s11[0],s12[0]};
Physical Line(10)={s21[0],s22[0]};
Physical Line(9)={-s12[3],-s22[3]};
Physical Line(8)={-s21[3],-s11[3]};


Physical Surface(20)={s11[1],s12[1],s21[1],s22[1]};
