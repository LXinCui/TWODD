input.E=30e9; % modulus
input.PR=0.25; % Poisson's ratio
input.PXX=0;
input.PYY=0;
input.PXY=0;   % Far-field stresses
%% 
clear rup;
stepSize=20;

for i=1:stepSize:4000  % rupture length
	if mod(i,100)==1;fprintf('Current step = %d\r\n',i);end
	N=100; % Number of elements
	input.X=[0.0:1/N:1.0]*i;
	input.Y=zeros(1,N+1);
	input.KOD(1:N)=3;   % KOD = 1 stress boundary conditions (sigmaN & sigmaS) in both normal and shear directions; 
						% KOD = 2 displacement discontinuity conditions (DN & DS) in both normal and shear directions; 
						% KOD = 3 DS/sigmaN in normal/shear direction;
						% KOD = 4 sigmaN/DS Normal/shear direction;						
	input.BVS(1:N)=1.0;
	input.BVN(1:N)=0.0;
	%input.FD.X;        % Calculate the stress and displacement at field points
	%input.FD.Y;
	clear BDoutput; clear FDoutput;
	[BDoutput,FDoutput]=TWODD_new(input); %BDoutput: outputs at boundary elements; FDoutput: outputs at field points
	rup.length(int16(i/stepSize)+1)=i;
	rup.strainEnergy(int16(i/stepSize)+1)=dot(diff(input.X),abs(BDoutput.SIGS)/2.0);
end 
figure;
plot(rup.length,rup.strainEnergy);
