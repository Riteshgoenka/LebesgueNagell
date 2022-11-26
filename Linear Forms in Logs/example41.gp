read("~/Downloads/lfl3-kit-main/lfl-utils-alpha1Variable.gp");

\\ L= 260, m=  5.0000, rho=  5.0000, chi= 7.5000, K= 2868204.725*logX, nonDegen log|Lambda|>-1.200211 e9*logX, nonDegenNUB=2.400423 e9, degenNUB1=    0.e-19, degenNUB2=2.411431 e9, degenNUB3=2.925711 e9, nUB=2.411431 e9, transB-b1

\\ L= 156, m=  5.0000, rho=  6.0000, chi=10.0000, K= 2530411.082*logX, nonDegen log|Lambda|>-7.072865 e8*logX, nonDegenNUB=1.414573 e9, degenNUB1=    0.e-19, degenNUB2=1.425139 e9, degenNUB3=1.753760 e9, nUB=1.425139 e9, transB-b1

\\ changed LB to 1.2*10^9
\\ L= 152, m=  6.5000, rho=  5.5000, chi=10.5000, K= 2647392.018*logX, nonDegen log|Lambda|>-6.859967 e8*logX, nonDegenNUB=1.371993 e9, degenNUB1=    0.e-19, degenNUB2=1.370891 e9, degenNUB3=1.691629 e9, nUB=1.371993 e9, transB-b1

\\ changed LB to 1.3*10^9
\\ L= 152, m=  6.7600, rho=  5.4000, chi=10.5600, K= 2643058.805*logX, nonDegen log|Lambda|>-6.775022 e8*logX, nonDegenNUB=1.355004 e9, degenNUB1=    0.e-19, degenNUB2=1.355367 e9, degenNUB3=1.650333 e9, nUB=1.355367 e9, transB-b1

eg41_search_it1(dbg=0) = {
	my(chiLB,chiUB,lLB,lUB,mLB,mUB,rhoLB,rhoUB,mStep,rhoStep,chiStep);

	lLB 	= 50;
	lUB 	= 350;
	mLB 	= 1; 	\\ m \geq 1
	mUB 	= 21;
	rhoLB 	= 2; 	\\ rho \geq 2
	rhoUB 	= 22;
	chiLB 	= 2;
	chiUB 	= 12;
	mStep 	= 1;
	rhoStep = 1;
	chiStep = 0.5;

	eg41_search_general(lLB,lUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB,mStep,rhoStep,chiStep,,dbg); 	\\ ",,dbg" here, as we are not passing in nUBInit (it defaults to 0)
}

eg41_search_it2(dbg=0) = {
	my(chiLB,chiUB,lLB,lUB,mLB,mUB,rhoLB,rhoUB,mStep,rhoStep,chiStep,nUBInit);

	lLB 	= 152;
	lUB 	= 152;
	mLB 	= 6; 	\\ m \geq 1
	mUB 	= 7;
	rhoLB 	= 5; 	\\ rho \geq 2
	rhoUB 	= 6;
	chiLB 	= 10;
	chiUB 	= 11;
	mStep 	= 0.04;
	rhoStep = 0.04;
	chiStep = 0.04;
	nUBInit = 1.371993*10^9;

	eg41_search_general(lLB,lUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB,mStep,rhoStep,chiStep,nUBInit,dbg); 	\\ ",,dbg" here, as we are not passing in nUBInit (it defaults to 0)
}

\\ L= 166, m=  6.0400, rho=  5.3600, chi=10.0200, K= 2556986.651*logX, nonDegen log|Lambda|>-7.126527 e8*logX, nonDegenNUB=1.425305 e9, degenNUB1=    0.e-19, degenNUB2=1.427511 e9, degenNUB3=1.732095 e9, nUB=1.427511 e9, transB-b1

eg41_check_it2(dbg=1) = {
	my(chi,bigL,m,rho,nUBInit);
	
	bigL 	= 166;
	m 		= 6.04;
	rho 	= 5.36;
	chi 	= 10.02;
	nUBInit = 2.600214*10^9;

	eg41_search_general(bigL,bigL,m,m,rho,rho,chi,chi,,,,nUBInit,dbg); 			\\ ",,dbg" here, as we are not passing in nUBInit (it defaults to 0)
}

eg41_search_general(bigLLB,bigLUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB,mStep=0.000001,rhoStep=0.000001,chiStep=0.000001,nUBInit=0,dbg=0) = {
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigD,bigK,d,d1,d2,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logW,logXLB,matveevChi,minNUB,nDegenUB,nLB,nNonDegenUB,nUB,val,w);

	bigD		= 2; 					\\ bigD=[Q(al_1,al_2):Q] -- used for Matveev's bounds
	matveevChi 	= 1; 					\\ matveevChi=[R(al_1,al_2):R] -- used for Matveev's bounds
	d 			= bigD/matveevChi; 		\\ d=[Q(al_1,al_2):Q]/[R(al_1,al_2):R] is the "degree" value used in the kit

	b1 			= n;
	b2 			= 2;
	b3 			= 2*(n+1); 				\\ since b3=2*r with 2<=r<=p+1
	nLB 		= 10^9;
	logXLB 		= 2*log(sqrt(nLB)-1);

	al3 		= 32 + 5*sqrt(41);
	hgtA3 		= (log(al3))/2;
	absLogA3	= abs(log(al3));

	al2 		= (7+sqrt(41))/(7-sqrt(41));
	hgtA2 		= (log(2)+log(al2))/2;
	absLogA2 	= abs(log(al2));

	al1 		= x; 					\\ this is just a placeholder so that al1=alpha_1 is considered to be a polynomial
	absLogA1 	= 2*log(al3) + (2*log(al3) + log(((sqrt(nLB)-1)^nLB+sqrt(41))/((sqrt(nLB)-1)^nLB-sqrt(41))) - 2*log(al2))/nLB;
	hgtA1 		= (logX+absLogA1)/2; 	\\ this must be correct though and consistent with logX usages that follow
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1 		= -logX/2;
	lamUB0 		= log(2.2*sqrt(41));

	if(nUBInit==0,
		nUBInit=get_matveev_ubnd(bigD,matveevChi,al1,absLogA1,hgtA1,al2,absLogA2,hgtA2,al3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,lamUB1,lamUB0,1);
	);
	minNUB=nUBInit;
	printf("used nUBInit=%9.6e\n",nUBInit);

	areBoundsOK=check_bounds(bigLLB,bigLUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB);
	if(areBoundsOK==0,
		return();
	);

	for(bigL=bigLLB,bigLUB, 	\\ L=5 is the lower bound in Theorem 4.1
	if(bigL%10==0,print("L=",bigL));
	forstep(m=mLB,mUB,mStep, 	\\ m \geq 1
		forstep(rho=rhoLB,rhoUB,rhoStep, 	\\ rho \geq 2
			a1 	= (rho-1)*absLogA1 + 2*bigD*hgtA1;
			a2 	= (rho-1)*absLogA2 + 2*bigD*hgtA2;
			a3 	= (rho-1)*absLogA3 + 2*bigD*hgtA3;
			if(dbg!=0,
				print("a1=",a1);
				print("a3=",a2);
				print("a3=",a3);
			);
			forstep(chi=chiLB,chiUB,chiStep,
				val=alpha1_check_params(d,al1,a1,absLogA1,hgtA1,al2,a2,absLogA2,hgtA2,al3,a3,absLogA3,hgtA3,b1,b2,b3,logXLB,nLB,bigL,m,rho,chi,nUBInit,lamUB1,lamUB0,dbg);
				minNUB=update_minNUB(val,bigL,m,rho,chi,minNUB,dbg);
			);
		);
	);
	);
}
