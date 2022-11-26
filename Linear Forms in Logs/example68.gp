read("~/Downloads/lfl3-kit-main/lfl-utils-alpha1Variable.gp");

\\ L= 210, m=  4.0000, rho=  6.0000, chi= 6.0000, K= 1235242.846*logX, nonDegen log|Lambda|>-4.647842 e8*logX, nonDegenNUB=9.295684 e8, degenNUB1=    0.e-19, degenNUB2=2.172008 e9, degenNUB3=9.481739 e8, nUB=9.481739 e8, transB-b1

\\ L= 142, m=  5.0000, rho=  6.0000, chi= 8.2500, K= 1044074.311*logX, nonDegen log|Lambda|>-2.656437 e8*logX, nonDegenNUB=5.312873 e8, degenNUB1=    0.e-19, degenNUB2=1.215315 e9, degenNUB3=5.340100 e8, nUB=5.340100 e8, transB-b1

\\ changed LB to 4.5*10^8
\\ L= 133, m=  5.5000, rho=  6.0000, chi= 8.5000, K= 1066299.999*logX, nonDegen log|Lambda|>-2.541036 e8*logX, nonDegenNUB=5.082071 e8, degenNUB1=    0.e-19, degenNUB2=1.158898 e9, degenNUB3=5.077975 e8, nUB=5.082071 e8, transB-b1

\\ changed LB to 4.9*10^8
\\ L= 133, m=  5.2500, rho=  6.1000, chi= 8.5500, K= 1050111.704*logX, nonDegen log|Lambda|>-2.525544 e8*logX, nonDegenNUB=5.051088 e8, degenNUB1=    0.e-19, degenNUB2=1.151468 e9, degenNUB3=5.045338 e8, nUB=5.051088 e8, transB-b1

eg68_search_it1(dbg=0) = {
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

	eg68_search_general(lLB,lUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB,mStep,rhoStep,chiStep,,dbg); 	\\ ",,dbg" here, as we are not passing in nUBInit (it defaults to 0)
}

eg68_check_it1(dbg=1) = {
	my(chi,bigL,m,rho);
	
	bigL 	= 233;
	m 		= 6.0;
	rho 	= 5.0;
	chi 	= 7.5;

	eg68_search_general(bigL,bigL,m,m,rho,rho,chi,chi,,,,,dbg); 			\\ ",,dbg" here, as we are not passing in nUBInit (it defaults to 0)
}

eg68_search_it2(dbg=0) = {
	my(chiLB,chiUB,lLB,lUB,mLB,mUB,rhoLB,rhoUB,mStep,rhoStep,chiStep,nUBInit);

	lLB 	= 133;
	lUB 	= 133;
	mLB 	= 5; 	\\ m \geq 1
	mUB 	= 6;
	rhoLB 	= 5; 	\\ rho \geq 2
	rhoUB 	= 7;
	chiLB 	= 8;
	chiUB 	= 9;
	mStep 	= 0.05;
	rhoStep = 0.1;
	chiStep = 0.05;
	nUBInit = 5.082071*10^8;

	eg68_search_general(lLB,lUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB,mStep,rhoStep,chiStep,nUBInit,dbg); 	\\ ",,dbg" here, as we are not passing in nUBInit (it defaults to 0)
}

eg68_check_it2(dbg=1) = {
	my(chi,bigL,m,rho,nUBInit);
	
	bigL 	= 49;
	m 		= 25.0;
	rho 	= 20.0;
	chi 	= 0.2;
	nUBInit = 1.883384*10^9;
	
	eg68_search_general(bigL,bigL,m,m,rho,rho,chi,chi,,,,nUBInit,dbg); 			\\ ",,dbg" here, as we are not passing in nUBInit (it defaults to 0)
}

eg68_search_general(bigLLB,bigLUB,mLB,mUB,rhoLB,rhoUB,chiLB,chiUB,mStep=0.000001,rhoStep=0.000001,chiStep=0.000001,nUBInit=0,dbg=0) = {
	my(a1,a2,a3,absLogA1,absLogA2,absLogA3,b1,b2,b3,bigD,bigK,d,d1,d2,hgtA1,hgtA2,hgtA3,lamUB0,lamUB1,logW,logXLB,matveevChi,minNUB,nDegenUB,nLB,nNonDegenUB,nUB,val,w);

	bigD		= 2; 					\\ bigD=[Q(al_1,al_2):Q] -- used for Matveev's bounds
	matveevChi 	= 1; 					\\ matveevChi=[R(al_1,al_2):R] -- used for Matveev's bounds
	d 			= bigD/matveevChi; 		\\ d=[Q(al_1,al_2):Q]/[R(al_1,al_2):R] is the "degree" value used in the kit

	b1 			= n;
	b2 			= 2;
	b3 			= 2*(n+1); 				\\ since b3=2*r with 2<=r<=p+1
	nLB 		= 3*10^8;
	logXLB 		= 2*log(sqrt(nLB)-1);

	al3 		= 4 + sqrt(17);
	hgtA3 		= (log(al3))/2;
	absLogA3	= abs(log(al3));

	al2 		= ((sqrt(17)+3)/(sqrt(17)-3))^2;
	hgtA2 		= (log(4)+log(al2))/2;
	absLogA2 	= abs(log(al2));

	al1 		= x; 					\\ this is just a placeholder so that al1=alpha_1 is considered to be a polynomial
	absLogA1 	= 2*log(al3);
	hgtA1 		= (logX+absLogA1)/2; 	\\ this must be correct though and consistent with logX usages that follow
	
	\\ assume that we have log |\Lambda| <lamUB1*n+lamUB0
	lamUB1 		= -logX/2;
	lamUB0 		= log(2.2*sqrt(68));

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
