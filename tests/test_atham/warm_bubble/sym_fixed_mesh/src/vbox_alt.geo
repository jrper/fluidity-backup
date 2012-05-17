Point(1)={500,1500,0};
Point(2)={500,0,0};

Point(3)={0   ,   0,0};
Point(4)={0   ,1500,0};
Point(5)={1000,1500,0};
Point(6)={1000,   0,0};



//N1=60;
//N3=10;
//N2=20;

N1=5;
N2=2;
N3=2;



//num[]=Extrude {0,-1240,0} {Point{1};Layers{N1};};
//num2[]=Extrude {0,260,0} {Point{2};Layers{N3};};
//s11[]=Extrude {500,0,0} {Line{num[1]};Layers{N2};};
//s12[]=Extrude {500,0,0} {Line{num2[1]};Layers{N2};};
//s21[]=Extrude {-500,0,0} {Line{num[1]};Layers{N2};};
//s22[]=Extrude {-500,0,0} {Line{num2[1]};Layers{N2};};

//num[]=Extrude {0,-1240,0} {Point{1};Layers{N1};Recombine;};
//num2[]=Extrude {0,260,0} {Point{2};Layers{N3};Recombine;};
//s11[]=Extrude {500,0,0} {Line{num[1]};Layers{N2};Recombine;};
//s12[]=Extrude {500,0,0} {Line{num2[1]};Layers{N2};Recombine;};
//s21[]=Extrude {-500,0,0} {Line{num[1]};Layers{N2};Recombine;};
//s22[]=Extrude {-500,0,0} {Line{num2[1]};Layers{N2};Recombine;};



//Physical Line(11)={s11[0],s12[0]};
//Physical Line(10)={s21[0],s22[0]};
//Physical Line(9)={-s12[3],-s22[3]};
//Physical Line(8)={-s21[3],-s11[3]};

Line(100)={3,4};
Line(101)={4,5};
Line(102)={5,6};
Line(103)={6,3};

Line Loop(200)={100,101,102,103};

Physical Line(8)={101};
Physical Line(9)={103};
Physical Line(10)={100};
Physical Line(11)={102};



Plane Surface(300)= {200};
//Physical Surface(20)={s11[1],s12[1],s21[1],s22[1]};


Transfinite Line{100,102}=90;
Transfinite Line{101,103}=60;
Transfinite Surface{300};

Recombine Surface{300};
Physical Surface(301)={300};