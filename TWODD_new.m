% TWODD.m - Translated from Fortran to MATLAB
% Main program for boundary element analysis
function [BDoutput,FDoutput]=TWODD_new(input)

global PI PR PR1 PR2 CON CONS
% S2
global SXXS SXXN SYYS SYYN SXYS SXYN UXS UXN UYS UYN
% S3
global C B D

% Initialize arrays
BDoutput=[];
FDoutput=[];

NUMBE=length(input.X)-1;
C = zeros(2*NUMBE,2*NUMBE);
B = zeros(2*NUMBE,1);
D = zeros(2*NUMBE,1);
COSBET = zeros(NUMBE,1);
SINBET = zeros(NUMBE,1);

E=input.E;
PR=input.PR;
KOD=input.KOD;
PYY=input.PYY;
PXX=input.PXX;
PXY=input.PXY;
XD=diff(input.X);
YD=diff(input.Y);
XM=(input.X(2:end)+input.X(1:end-1))./2;
YM=(input.Y(2:end)+input.Y(1:end-1))./2;

SW=sqrt(XD.^2 + YD.^2);
A=0.5*SW;
SINBET=YD./SW;
COSBET=XD./SW;

PI = 4*atan(1.0);
CON = 1/(4*PI*(1-PR));
CONS = E/(1+PR);
PR1 = 1-2*PR;
PR2 = 2*(1-PR);

for i=1:NUMBE
	MN = 2*i;
	MS = MN-1;
	B(MS) = input.BVS(i);
	B(MN) = input.BVN(i);
end

% ADJUST STRESS BOUNDARY VALUES TO ACCOUNT FOR INITIAL STRESSES
for N = 1:NUMBE
    NN = 2*N;
    NS = NN-1;
    COSB = COSBET(N);
    SINB = SINBET(N);
    SIGS = (PYY-PXX)*SINB*COSB + PXY*(COSB^2 - SINB^2);
    SIGN = PXX*SINB^2 - 2*PXY*SINB*COSB + PYY*COSB^2;
    switch KOD(N)
        case 1
            B(NS) = B(NS) - SIGS;
            B(NN) = B(NN) - SIGN;
        case 2
            B(NN) = B(NN) - SIGN;
        case 3
            B(NS) = B(NS) - SIGS;
    end
end

% COMPUTE INFLUENCE COEFFICIENTS AND SET UP SYSTEM OF ALGEBRAIC EQUATIONS
KSYM=1;
for I = 1:NUMBE
    IN = 2*I;
    IS = IN-1;
    XI = XM(I);
    YI = YM(I);
    COSBI = COSBET(I);
    SINBI = SINBET(I);
    KODE = KOD(I);
    for J = 1:NUMBE
        JN = 2*J;
        JS = JN-1;
        INITL();
        XJ = XM(J);
        YJ = YM(J);
        COSBJ = COSBET(J);
        SINBJ = SINBET(J);
        AJ = A(J);
        switch KSYM
            case 1
                COEFF(XI,YI,XJ,YJ,AJ,COSBJ,SINBJ,1);
            case 2
                XJ2 = 2*XSYM - XM(J);
                COEFF(XI,YI,XJ2,YJ,AJ,COSBJ,-SINBJ,-1);
            case 3
                YJ2 = 2*YSYM - YM(J);
                COEFF(XI,YI,XJ,YJ2,AJ,-COSBJ,SINBJ,-1);
            case 4
                XJ2 = 2*XSYM - XM(J);
                COEFF(XI,YI,XJ2,YJ,AJ,COSBJ,-SINBJ,-1);
                COEFF(XI,YI,XJ,YJ2,AJ,-COSBJ,SINBJ,-1);
                COEFF(XI,YI,XJ2,YJ2,AJ,-COSBJ,-SINBJ,1);
        end
        switch KODE
            case 1
                C(IS,JS) = (SYYS-SXXS)*SINBI*COSBI + SXYS*(COSBI^2-SINBI^2);
                C(IS,JN) = (SYYN-SXXN)*SINBI*COSBI + SXYN*(COSBI^2-SINBI^2);
                C(IN,JS) = SXXS*SINBI^2 - 2*SXYS*SINBI*COSBI + SYYS*COSBI^2;
                C(IN,JN) = SXXN*SINBI^2 - 2*SXYN*SINBI*COSBI + SYYN*COSBI^2;
            case 2
                C(IS,JS) = UXS*COSBI + UYS*SINBI;
                C(IS,JN) = UXN*COSBI + UYN*SINBI;
                C(IN,JS) = -UXS*SINBI + UYS*COSBI;
                C(IN,JN) = -UXN*SINBI + UYN*COSBI;
            case 3
                C(IS,JS) = UXS*COSBI + UYS*SINBI;
                C(IS,JN) = UXN*COSBI + UYN*SINBI;
                C(IN,JS) = SXXS*SINBI^2 - 2*SXYS*SINBI*COSBI + SYYS*COSBI^2;
                C(IN,JN) = SXXN*SINBI^2 - 2*SXYN*SINBI*COSBI + SYYN*COSBI^2;
            case 4
                C(IS,JS) = (SYYS-SXXS)*SINBI*COSBI + SXYS*(COSBI^2-SINBI^2);
                C(IS,JN) = (SYYN-SXXN)*SINBI*COSBI + SXYN*(COSBI^2-SINBI^2);
                C(IN,JS) = -UXS*SINBI + UYS*COSBI;
                C(IN,JN) = -UXN*SINBI + UYN*COSBI;
        end
    end
end

% SOLVE SYSTEM OF ALGEBRAIC EQUATIONS
%N = 2*NUMBE;
%D = SOLVE(C, B, N);
D = C\B;

for i=1:NUMBE
	BDoutput.DN(i)=D(2*i);
	BDoutput.DS(i)=D(2*i-1);
end
BDoutput.DN=BDoutput.DN';
BDoutput.DS=BDoutput.DS';
	

% COMPUTE BOUNDARY DISPLACEMENTS AND STRESSES
%fprintf(fid_out, '\n     DISPLACEMENTS AND STRESSES AT MIDPOINTS OF BOUNDARY ELEMENTS.\n\n   ELEMENT        DS     US(-)     US(+)     DN     UN(-)     UN(+)     UX(-)     UY(-)     UX(+)     UY(+)   SIGMA-S   SIGMA-N\n');
for I = 1:NUMBE
    IN = 2*I;
    IS = IN-1;
    XI = XM(I);
    YI = YM(I);
    COSBI = COSBET(I);
    SINBI = SINBET(I);
    UXNEG = 0; UYNEG = 0; SIGXX = PXX; SIGYY = PYY; SIGXY = PXY;
    for J = 1:NUMBE
        JN = 2*J;
        JS = JN-1;
        INITL();
        XJ = XM(J);
        YJ = YM(J);
        AJ = A(J);
        COSBJ = COSBET(J);
        SINBJ = SINBET(J);
        switch KSYM
            case 1
                COEFF(XI,YI,XJ,YJ,AJ,COSBJ,SINBJ,1);
            case 2
                XJ2 = 2*XSYM - XM(J);
                COEFF(XI,YI,XJ2,YJ,AJ,COSBJ,-SINBJ,-1);
            case 3
                YJ2 = 2*YSYM - YM(J);
                COEFF(XI,YI,XJ,YJ2,AJ,-COSBJ,SINBJ,-1);
            case 4
                XJ2 = 2*XSYM - XM(J);
                COEFF(XI,YI,XJ2,YJ,AJ,COSBJ,-SINBJ,-1);
                COEFF(XI,YI,XJ,YJ2,AJ,-COSBJ,SINBJ,-1);
                COEFF(XI,YI,XJ2,YJ2,AJ,-COSBJ,-SINBJ,1);
        end
        UXNEG = UXNEG + UXS*D(JS) + UXN*D(JN);
        UYNEG = UYNEG + UYS*D(JS) + UYN*D(JN);
        SIGXX = SIGXX + SXXS*D(JS) + SXXN*D(JN);
        SIGYY = SIGYY + SYYS*D(JS) + SYYN*D(JN);
        SIGXY = SIGXY + SXYS*D(JS) + SXYN*D(JN);
    end
    USNEG = UXNEG*COSBI + UYNEG*SINBI;
    UNNEG = -UXNEG*SINBI + UYNEG*COSBI;
    USPOS = USNEG - D(IS);
    UNPOS = UNNEG - D(IN);
    UXPOS = USPOS*COSBI - UNPOS*SINBI;
    UYPOS = USPOS*SINBI + UNPOS*COSBI;
    SIGS = (SIGYY-SIGXX)*SINBI*COSBI + SIGXY*(COSBI^2-SINBI^2);
    SIGN = SIGXX*SINBI^2 - 2*SIGXY*SINBI*COSBI + SIGYY*COSBI^2;
    %fprintf(fid_out, '%10d%15.3e%15.3e%15.3e%15.3e%15.3e%15.3e%15.3e%15.3e%15.3e%15.3e%15.3e%15.3e\n', I, D(IS), USNEG, USPOS, D(IN), UNNEG, UNPOS, UXNEG, UYNEG, UXPOS, UYPOS, SIGS, SIGN);
	BDoutput.SIGS(I)=SIGS;
	BDoutput.SIGN(I)=SIGN;
end
BDoutput.SIGS=BDoutput.SIGS';
BDoutput.SIGN=BDoutput.SIGN';

% COMPUTE DISPLACEMENTS AND STRESSES AT SPECIFIED POINTS IN BODY
if ~isfield('FD',input) 
	FDoutput=[];
else 
    %fprintf(fid_out, '\n     DISPLACEMENTS AND STRESSES AT SPECIFIED POINTS IN THE BODY.\n\n    POINT    X CO-ORD    Y CO-ORD          UX          UY      SIGXX      SIGYY      SIGXY\n');
    NPOINT = length(input.FD.X);
    for N = 1:NPOINT
		XP=input.FD.X(N);
		YP=input.FD.Y(N);
		UX = 0; UY = 0; SIGXX = PXX; SIGYY = PYY; SIGXY = PXY;
		for J = 1:NUMBE
			JN = 2*J;
			JS = JN-1;
			INITL();
			XJ = XM(J);
			YJ = YM(J);
			AJ = A(J);
			if sqrt((XP-XJ)^2 + (YP-YJ)^2) < 2*AJ
				continue;
			end
			COSBJ = COSBET(J);
			SINBJ = SINBET(J);
			switch KSYM
				case 1
					COEFF(XP,YP,XJ,YJ,AJ,COSBJ,SINBJ,1);
				case 2
					XJ2 = 2*XSYM - XM(J);
					COEFF(XP,YP,XJ2,YJ,AJ,COSBJ,-SINBJ,-1);
				case 3
					YJ2 = 2*YSYM - YM(J);
					COEFF(XP,YP,XJ,YJ2,AJ,-COSBJ,SINBJ,-1);
				case 4
					XJ2 = 2*XSYM - XM(J);
					COEFF(XP,YP,XJ2,YJ,AJ,COSBJ,-SINBJ,-1);
					COEFF(XP,YP,XJ,YJ2,AJ,-COSBJ,SINBJ,-1);
					COEFF(XP,YP,XJ2,YJ2,AJ,-COSBJ,-SINBJ,1);
			end
			UX = UX + UXS*D(JS) + UXN*D(JN);
			UY = UY + UYS*D(JS) + UYN*D(JN);
			SIGXX = SIGXX + SXXS*D(JS) + SXXN*D(JN);
			SIGYY = SIGYY + SYYS*D(JS) + SYYN*D(JN);
			SIGXY = SIGXY + SXYS*D(JS) + SXYN*D(JN);
		end
		FDoutput.UX(N)=UX;
		FDoutput.UY(N)=UY;
		FDoutput.SIGXX(N)=SIGXX;
		FDoutput.SIGYY(N)=SIGYY;
		FDoutput.SIGXY(N)=SIGXY;
		%fprintf(fid_out, '%9d%13.3e%13.3e%13.3e%13.3e%13.3e%13.3e%13.3e\n', NPOINT, XP, YP, UX, UY, SIGXX, SIGYY, SIGXY);
    end
end
end 
% --- Subroutines as MATLAB functions ---

function INITL()
    global SXXS SXXN SYYS SYYN SXYS SXYN UXS UXN UYS UYN
    SXXS=0; SXXN=0; SYYS=0; SYYN=0; SXYS=0; SXYN=0;
    UXS=0; UXN=0; UYS=0; UYN=0;
end

function COEFF(X,Y,CX,CY,A,COSB,SINB,MSYM)
    global PI PR PR1 PR2 CON CONS
    global SXXS SXXN SYYS SYYN SXYS SXYN UXS UXN UYS UYN
    COS2B = COSB^2 - SINB^2;
    SIN2B = 2*SINB*COSB;
    COSB2 = COSB^2;
    SINB2 = SINB^2;
    XB = (X-CX)*COSB + (Y-CY)*SINB;
    YB = -(X-CX)*SINB + (Y-CY)*COSB;
    R1S = (XB-A)^2 + YB^2;
    R2S = (XB+A)^2 + YB^2;
    FL1 = 0.5*log(R1S);
    FL2 = 0.5*log(R2S);
    FB2 = CON*(FL1-FL2);
    if YB ~= 0
        FB3 = -CON*(atan((XB+A)/YB) - atan((XB-A)/YB));
    else
        FB3 = 0;
        if abs(XB) < A
            FB3 = CON*PI;
        end
    end
    FB4 = CON*(YB/R1S - YB/R2S);
    FB5 = CON*((XB-A)/R1S - (XB+A)/R2S);
    FB6 = CON*(((XB-A)^2-YB^2)/R1S^2 - ((XB+A)^2-YB^2)/R2S^2);
    FB7 = 2*CON*YB*((XB-A)/R1S^2 - (XB+A)/R2S^2);
    UXDS = -PR1*SINB*FB2 + PR2*COSB*FB3 + YB*(SINB*FB4 - COSB*FB5);
    UXDN = -PR1*COSB*FB2 - PR2*SINB*FB3 - YB*(COSB*FB4 + SINB*FB5);
    UYDS = PR1*COSB*FB2 + PR2*SINB*FB3 - YB*(COSB*FB4 + SINB*FB5);
    UYDN = -PR1*SINB*FB2 + PR2*COSB*FB3 - YB*(SINB*FB4 - COSB*FB5);
    SXXDS = CONS*(2*COSB2*FB4 + SIN2B*FB5 + YB*(COS2B*FB6 - SIN2B*FB7));
    SXXDN = CONS*(-FB5 + YB*(SIN2B*FB6 + COS2B*FB7));
    SYYDS = CONS*(2*SINB2*FB4 - SIN2B*FB5 - YB*(COS2B*FB6 - SIN2B*FB7));
    SYYDN = CONS*(-FB5 - YB*(SIN2B*FB6 + COS2B*FB7));
    SXYDS = CONS*(SIN2B*FB4 - COS2B*FB5 + YB*(SIN2B*FB6 + COS2B*FB7));
    SXYDN = CONS*(-YB*(COS2B*FB6 - SIN2B*FB7));
    UXS = UXS + MSYM*UXDS;
    UXN = UXN + UXDN;
    UYS = UYS + MSYM*UYDS;
    UYN = UYN + UYDN;
    SXXS = SXXS + MSYM*SXXDS;
    SXXN = SXXN + SXXDN;
    SYYS = SYYS + MSYM*SYYDS;
    SYYN = SYYN + SYYDN;
    SXYS = SXYS + MSYM*SXYDS;
    SXYN = SXYN + SXYDN;
end

function D = SOLVE(A, B, N)
    % Gaussian elimination with back substitution
    A = A(1:N,1:N);
    B = B(1:N);
    for J = 1:N-1
        L = J+1;
        for JJ = L:N
            XM = A(JJ,J)/A(J,J);
            for I = J:N
                A(JJ,I) = A(JJ,I) - A(J,I)*XM;
            end
            B(JJ) = B(JJ) - B(J)*XM;
        end
    end
    D = zeros(N,1);
    D(N) = B(N)/A(N,N);
    for J = 1:N-1
        JJ = N-J;
        L = JJ+1;
        SUM = 0;
        for I = L:N
            SUM = SUM + A(JJ,I)*D(I);
        end
        D(JJ) = (B(JJ) - SUM)/A(JJ,JJ);
    end
end
