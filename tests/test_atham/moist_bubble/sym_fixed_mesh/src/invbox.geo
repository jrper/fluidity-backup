Point(1)={500,1500,0};

N1=60;
N2=20;


num[]=Extrude {0,-1500,0} {Point{1};Layers{N1};};
s1[]=Extrude {500,0,0} {Line{num[1]};Layers{N2};};
s2[]=Extrude {-500,0,0} {Line{num[1]};Layers{N2};};

Physical Line(11)={s1[0]};
Physical Line(10)={s2[0]};
Physical Line(9)={s1[2],-s2[2]};
Physical Line(8)={s1[3],-s2[3]};


Physical Surface(20)={s1[1],s2[1]};
