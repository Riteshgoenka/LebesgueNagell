read("~/Downloads/lfl3-kit-main/lfl-utils-alpha1Variable.gp");

\\ L= 247, m=  6.0000, rho=  5.0000, chi= 8.5000, K=10222492.143*logX, nonDegen log|Lambda|>-4.063759 e9*logX, nonDegenNUB=8.127518 e9, degenNUB1=    0.e-19, degenNUB2=8.072696 e9, degenNUB3=1.537479 e10, nUB=8.127518 e9, transB-b1

\\ L= 180, m=  7.0000, rho=  5.0000, chi=11.0000, K= 8691187.652*logX, nonDegen log|Lambda|>-2.517827 e9*logX, nonDegenNUB=5.035654 e9, degenNUB1=    0.e-19, degenNUB2=4.946913 e9, degenNUB3=9.506115 e9, nUB=5.035654 e9, transB-b1

\\ L= 180, m=  6.6000, rho=  5.0400, chi=10.8800, K= 8339051.040*logX, nonDegen log|Lambda|>-2.427774 e9*logX, nonDegenNUB=4.855547 e9, degenNUB1=    0.e-19, degenNUB2=4.846831 e9, degenNUB3=9.358625 e9, nUB=4.855547 e9, transB-b1

\\ changed LB to 4.5*10^9
\\ L= 180, m=  6.3200, rho=  5.1200, chi=10.8800, K= 7878790.777*logX, nonDegen log|Lambda|>-2.316111 e9*logX, nonDegenNUB=4.632222 e9, degenNUB1=    0.e-19, degenNUB2=4.638879 e9, degenNUB3=8.843229 e9, nUB=4.638879 e9, transB-b1

eg73_search_it1(dbg=0) = {
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

	eg73_search_general(lLB,lUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB,mStep,rhoStep,chiStep,,dbg); 	\\ ",,dbg" here, as we are not passing in nUBInit (it defaults to 0)
}

eg73_search_it2(dbg=0) = {
	my(chiLB,chiUB,lLB,lUB,mLB,mUB,rhoLB,rhoUB,mStep,rhoStep,chiStep,nUBInit);

	lLB 	= 180;
	lUB 	= 180;
	mLB 	= 6; 	\\ m \geq 1
	mUB 	= 8;
	rhoLB 	= 4; 	\\ rho \geq 2
	rhoUB 	= 6;
	chiLB 	= 10;
	chiUB 	= 12;
	mStep 	= 0.04;
	rhoStep = 0.04;
	chiStep = 0.04;
	nUBInit = 4.855547*10^9;

	eg73_search_general(lLB,lUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB,mStep,rhoStep,chiStep,nUBInit,dbg); 	\\ ",,dbg" here, as we are not passing in nUBInit (it defaults to 0)
}

eg73_check_it1(dbg=1) = {
	my(chi,bigL,m,rho);
	
	bigL 	= 247;
	m 		= 6.0;
	rho 	= 5.0;
	chi 	= 8.5;

	eg73_search_general(bigL,bigL,m,m,rho,rho,chi,chi,,,,,dbg); 			\\ ",,dbg" here, as we are not passing in nUBInit (it defaults to 0)
}

eg73_check_it2(dbg=1) = {
	my(chi,bigL,m,rho,nUBInit);
	
	bigL 	= 180;
	m 		= 7.0;
	rho 	= 5.0;
	chi 	= 11.0;
	nUBInit = 8.127518*10^9;
	
	eg73_search_general(bigL,bigL,m,m,rho,rho,chi,chi,,,,nUBInit,dbg); 			\\ ",,dbg" here, as we are not passing in nUBInit (it defaults to 0)
}

eg73_check_it3(dbg=1) = {
	my(chi,bigL,m,rho,nUBInit);
	
	bigL 	= 180;
	m 		= 6.6;
	rho 	= 5.04;
	chi 	= 10.88;
	nUBInit = 5.035654*10^9;
	
	eg73_search_general(bigL,bigL,m,m,rho,rho,chi,chi,,,,nUBInit,dbg); 			\\ ",,dbg" here, as we are not passing in nUBInit (it defaults to 0)
}

eg73_check_it4(dbg=1) = {
	my(chi,bigL,m,rho,nUBInit);
	
	bigL 	= 180;
	m 		= 6.32;
	rho 	= 5.12;
	chi 	= 10.88;
	nUBInit = 4.855547*10^9;
	
	eg73_search_general(bigL,bigL,m,m,rho,rho,chi,chi,,,,nUBInit,dbg); 			\\ ",,dbg" here, as we are not passing in nUBInit (it defaults to 0)
}

eg73_search_general(bigLLB,bigLUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB,mStep=0.000001,rhoStep=0.000001,chiStep=0.000001,nUBInit=0,dbg=0) = {
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigD,bigK,d,d1,d2,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logW,logXLB,matveevChi,minNUB,nDegenUB,nLB,nNonDegenUB,nUB,val,w);

	bigD		= 2; 					\\ bigD=[Q(al_1,al_2):Q] -- used for Matveev's bounds
	matveevChi 	= 1; 					\\ matveevChi=[R(al_1,al_2):R] -- used for Matveev's bounds
	d 			= bigD/matveevChi; 		\\ d=[Q(al_1,al_2):Q]/[R(al_1,al_2):R] is the "degree" value used in the kit

	b1 			= n;
	b2 			= 2;
	b3 			= 2*(n+1); 				\\ since b3=2*r with 2<=r<=p+1
	nLB 		= 10^9;
	logXLB 		= 2*log(sqrt(nLB)-1);

	al3 		= 1068 + 125*sqrt(73);
	hgtA3 		= (log(al3))/2;
	absLogA3	= abs(log(al3));

	al2 		= (9+sqrt(73))/(9-sqrt(73));
	hgtA2 		= (log(2)+log(al2))/2;
	absLogA2 	= abs(log(al2));

	al1 		= x; 					\\ this is just a placeholder so that al1=alpha_1 is considered to be a polynomial
	absLogA1 	= 2*log(al3) + (2*log(al3) + log(((sqrt(nLB)-1)^nLB+sqrt(73))/((sqrt(nLB)-1)^nLB-sqrt(73))) - 2*log(al2))/nLB;
	hgtA1 		= (logX+absLogA1)/2; 	\\ this must be correct though and consistent with logX usages that follow
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1 		= -logX/2;
	lamUB0 		= log(2.2*sqrt(73));

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
