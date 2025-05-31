input.E=30e9;
input.PR=0.25;
input.PXX=0;
input.PYY=0;
input.PXY=0;
input.X=[0:1/50:1];
input.Y=zeros(1,51);
input.KOD(1:50)=2;
input.BVS(1:50)=1.0;
input.BVN(1:50)=0.0;
%input.FD.X;
%input.FD.Y;

[BDoutput,FDoutput]=TWODD_new(input);