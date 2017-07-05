#include "ab_initio_potential.h"
#include <math.h>
#include <stdio.h>

double * legendre_array(int N, double cosT)
{
	double * leg = (double*)malloc((N+1)*sizeof(double));
	leg[0] = 1;	
	leg[1] = cosT;

	for (int i =1; i < N; ++i)
	{
		leg[i+1] = ((2.0*i+1.0)/(i+1.0))*cosT*leg[i] - ((double)i/(i+1.0))*leg[i-1];
	}
	return leg;

}

double ab_initio_pot(double R, double Theta)
{
	double t2, t4, t5, t6, t10, t11, t14, t15, t16, t18, t20, t21, t22, t23, t24, t25, t26, t29, t31, t39, t47, t48;
	double cosT = cos(Theta);

	double *legP = legendre_array(10,cosT);

		if (R<4.500000)
	{
	t2 = exp(-0.15587635029119e1 * R);
	t5 = exp(-0.13910571666810e1 * R);
	t6 = legP[2];
	t10 = exp(-0.11314491915437e1 * R);
	t11 = legP[4];
	t15 = exp(-0.63429130737700e0 * R);
	t16 = legP[6];
	t20 = exp(0.49580080306684e-1 * R);
	t21 = legP[8];
	t25 = exp(-0.10290558301028e-1 * R);
	t26 = legP[10];
	t48 = 0.41801730892757e8 * t2 + 0.48642926965367e8 * t5 * t6 + 0.74492214522289e7 * t10 * t11 + 0.21822899669372e6 * t15 * t16 + 0.20109486940810e4 * 
			t20 * t21 + 0.60812928463236e3 * t25 * t26;

			
	}
	else if (R<4.750000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = 0.45717035866570e7 - 0.26840768764320e7 * R + 0.53421446026433e6 * t2 - 0.35924592941254e5 * t4 + 0.11057970695282e8 * t6 - 0.65871285639296e7 * 
			t6 * R + 0.13319170648237e7 * t6 * t2 - 0.91020264029875e5 * t6 * t4 + 0.69199195771701e7 * t14 - 0.42874360880077e7 * t14 * R + 0.89865578933964e6 *
			t14 * t2 - 0.63412172257467e5 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = 0.34466709030487e7 * t23 - 0.22234799340508e7 * t23 * R + 0.48122911483743e6 * t23 * t2 - 0.34823942580519e5 * t23 * t4 + 0.13997106233326e7 * 
			t31 - 0.91535265960639e6 * t31 * R + 0.19980355866819e6 * t31 * t2 - 0.14530942284715e5 * t31 * t4 + 0.23117962056184e6 * t39 - 0.15102384289339e6 * 
			t39 * R + 0.32960219368772e5 * t39 * t2 - 0.23971046293147e4 * t39 * t4;
	t48 = t22 + t47;
	}
	else if (R<5.000000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = 0.62793033600858e6 - 0.19327271812782e6 * R + 0.98346374634341e4 * t2 + 0.87399111494903e3 * t4 - 0.11719558775144e7 * t6 + 0.11370355873102e7 * t6 * R - 0.29422275648990e6 * t6 * t2 + 0.23094811150029e5 * t6 * t4 - 0.38995099806608e7 * t14 + 0.25458878432540e7 * t14 * R - 0.53993872250494e6 * t14 * t2 + 0.37541828573731e5 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = -0.30733947099906e7 * t23 + 0.18944562426055e7 * t23 * R - 0.38570481709023e6 * t23 * t2 + 0.26013526326683e5 * t23 * t4 - 0.11458946884890e7 * t31 + 0.69239806364938e6 * t31 * R - 0.13867027780672e6 * t31 * t2 + 0.92216076433483e4 * t31 * t4 - 0.11359508779163e6 * t39 + 0.66728604487747e5 * t39 * R - 0.12882401132519e5 * t39 * t2 + 0.81992137077595e3 * t39 * t4;
	t48 = t22 + t47;
	}
	else if (R<5.250000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = 0.19318024162470e7 - 0.97559596627091e6 * R + 0.16629928709205e6 * t2 - 0.95569855269585e4 * t4 + 0.41230152545009e7 * t6 - 0.20399470918990e7 * t6 * R + 0.34117377935192e6 * t6 * t2 - 0.19264957906092e5 * t6 * t4 + 0.15148265023583e7 * t14 - 0.70271404655758e6 * t14 * R + 0.10978165545737e6 * t14 * t2 - 0.57728632904225e4 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = 0.55595999271146e4 * t23 + 0.47083656654884e5 * t23 * R - 0.16230299900107e5 * t23 * t2 + 0.13818918473423e4 * t23 * t4 - 0.28771797467766e6 * t31 + 0.17749203536256e6 * t31 * R - 0.35689072149351e5 * t31 * t2 + 0.23561939328573e4 * t31 * t4 - 0.19370835069261e6 * t39 + 0.11479656222834e6 * t39 * R - 0.22495992680637e5 * t39 * t2 + 0.14608274739838e4 * t39 * t4;
	t48 = t22 + t47;
	}
	else if (R<5.500000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = 0.14195393368546e7 - 0.68287420661810e6 * R + 0.11054276144390e6 * t2 - 0.60168886604093e4 * t4 + 0.31557474284718e7 * t6 - 0.14872226198823e7 * t6 * R + 0.23589292753924e6 * t6 * t2 - 0.12580459378303e5 * t6 * t4 + 0.13351759877765e7 * t14 - 0.60005660965372e6 * t14 * R + 0.90227857951869e5 * t14 * t2 - 0.45313523376925e4 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = 0.25760893184280e6 * t23 - 0.96944533011222e5 * t23 * R + 0.11203640988675e5 * t23 * t2 - 0.35994566940579e3 * t23 * t4 + 0.49717421764701e5 * t31 - 0.15328191175934e5 * t31 * R + 0.10385900484567e4 * t31 * t2 + 0.24278872678960e2 * t31 * t4 + 0.47505158307385e5 * t39 - 0.23039728628803e5 * t39 * R + 0.37585389111991e4 * t39 * t2 - 0.20612691279945e3 * t39 * t4;
	t48 = t22 + t47;
	}
	else if (R<5.750000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = 0.12629721014668e7 - 0.59747389640657e6 * R + 0.95015432314530e5 * t2 - 0.50758384101445e4 * t4 + 0.30992754717798e7 * t6 - 0.14564197344140e7 * t6 * R + 0.23029240290863e6 * t6 * t2 - 0.12241033643114e5 * t6 * t4 + 0.16636615861008e7 * t14 - 0.77923057237607e6 * t14 * R + 0.12280494208321e6 * t14 * t2 - 0.65057210729250e4 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = 0.55259653517478e6 * t23 - 0.25784686210140e6 * t23 * R + 0.40458609914160e5 * t23 * t2 - 0.21329740891322e4 * t23 * t4 + 0.14466770974384e6 * t31 - 0.67119257346375e5 * t31 * R + 0.10455147533991e5 * t31 * t2 - 0.54642158098979e3 * t31 * t4 + 0.33772177588958e5 * t39 - 0.15549011873297e5 * t39 * R + 0.23965904101980e4 * t39 * t2 - 0.12358457940545e3 * t39 * t4;
	t48 = t22 + t47;
	}
	else if (R<6.000000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = 0.94978948289215e6 - 0.43407426932413e6 * R + 0.66598105865409e5 * t2 - 0.34284571667173e4 * t4 + 0.23855365192046e7 * t6 - 0.10840341939400e7 * t6 * R + 0.16552970021750e6 * t6 * t2 - 0.84866740668168e4 * t6 * t4 + 0.13105295736197e7 * t14 - 0.59498778325552e6 * t14 * R + 0.90762717888324e5 * t14 * t2 - 0.46482008297436e4 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = 0.45526562948209e6 * t23 - 0.20706552000086e6 * t23 * R + 0.31627072157544e5 * t23 * t2 - 0.16210008858502e4 * t23 * t4 + 0.13303883578822e6 * t31 - 0.61052018760835e5 * t31 * R + 0.93999756060713e4 * t31 * t2 - 0.48525219386401e3 * t31 * t4 + 0.37884633507102e5 * t39 - 0.17694641047982e5 * t39 * R + 0.27697433101431e4 * t39 * t2 - 0.14521663157618e3 * t39 * t4;
	t48 = t22 + t47;
	}
	else if (R<6.250000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = 0.71770416128376e6 - 0.31803160851991e6 * R + 0.47257662398039e5 * t2 - 0.23539880851967e4 * t4 + 0.18449542347862e7 * t6 - 0.81374305173082e6 * t6 * R + 0.12048117651597e6 * t6 * t2 - 0.59839783056209e4 * t6 * t4 + 0.10334995862982e7 * t14 - 0.45647278959472e6 * t14 * R + 0.67676885611527e5 * t14 * t2 - 0.33656545921438e4 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = 0.36507495437230e6 * t23 - 0.16197018244596e6 * t23 * R + 0.24111182565063e5 * t23 * t2 - 0.12034514640456e4 * t23 * t4 + 0.10340794082099e6 * t31 - 0.46236571277218e5 * t31 * R + 0.69307343588017e4 * t31 * t2 - 0.34807212457125e3 * t31 * t4 + 0.25117136477357e5 * t39 - 0.11310892533109e5 * t39 * R + 0.17057852243309e4 * t39 * t2 - 0.86107849031062e2 * t39 * t4;
	t48 = t22 + t47;
	}
	else if (R<6.500000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = 0.51947929515293e6 - 0.22288367277713e6 * R + 0.32033992679193e5 * t2 - 0.15420590335249e4 * t4 + 0.13532936743810e7 * t6 - 0.57774598273630e6 * t6 * R + 0.82721645476849e5 * t6 * t2 - 0.39701366502011e4 * t6 * t4 + 0.75351409297823e6 * t14 - 0.32207975280114e6 * t14 * R + 0.46173999724556e5 * t14 * t2 - 0.22188340115053e4 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = 0.26123082864649e6 * t23 - 0.11212500209757e6 * t23 * R + 0.16135953709321e5 * t23 * t2 - 0.77810592507268e3 * t23 * t4 + 0.72089505249899e5 * t31 - 0.31203722203094e5 * t31 * R + 0.45254785069419e4 * t31 * t2 - 0.21979181247207e3 * t31 * t4 + 0.17694017350378e5 * t39 - 0.77477953521589e4 * t39 * R + 0.11356896753790e4 * t39 * t2 - 0.55702753086958e2 * t39 * t4;
	t48 = t22 + t47;
	}
	else if (R<6.750000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = 0.37065694592956e6 - 0.15419643467403e6 * R + 0.21466725278717e5 * t2 - 0.10001478847826e4 * t4 + 0.97928432685695e6 * t6 - 0.40512628387905e6 * t6 * R + 0.56164768729581e5 * t6 * t2 - 0.26082455349566e4 * t6 * t4 + 0.54138985882135e6 * t14 - 0.22417626011337e6 * t14 * R + 0.31111923926434e5 * t14 * t2 - 0.14464198680119e4 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = 0.18320841990431e6 * t23 - 0.76114659601186e5 * t23 * R + 0.10595901017569e5 * t23 * t2 - 0.49400065882898e3 * t23 * t4 + 0.47427892385178e5 * t31 - 0.19821439342453e5 * t31 * R + 0.27743580668435e4 * t31 * t2 - 0.12999076426189e3 * t31 * t4 + 0.10367695988582e5 * t39 - 0.43664162620988e4 * t39 * R + 0.61547750767745e3 * t39 * t2 - 0.29025206025339e2 * t39 * t4;
	t48 = t22 + t47;
	}
	else if (R<7.000000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = 0.25884598362595e6 - 0.10450267365021e6 * R + 0.14104686608521e5 * t2 - 0.63659041958770e3 * t4 + 0.69419983092550e6 * t6 - 0.27842206346506e6 * t6 * R + 0.37393773112694e5 * t6 * t2 - 0.16812827884435e4 * t6 * t4 + 0.38043259681779e6 * t14 - 0.15263969922289e6 * t14 * R + 0.20513914905623e5 * t14 * t2 - 0.92306139784835e3 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = 0.12601282110069e6 * t23 - 0.50694393466241e5 * t23 * R + 0.68299356642436e4 * t23 * t2 - 0.30802706113391e3 * t23 * t4 + 0.31067033773137e5 * t31 - 0.12549946625991e5 * t31 * R + 0.16970998866269e4 * t31 * t2 - 0.76792829436374e2 * t31 * t4 + 0.62802993461365e4 * t39 - 0.25497955321232e4 * t39 * R + 0.34634851064403e3 * t39 * t2 - 0.15734885184182e2 * t39 * t4;
	t48 = t22 + t47;
	}
	else if (R<7.250000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = 0.17761237991913e6 - 0.69688272061566e5 * R + 0.91312006672867e4 * t2 - 0.39975775571942e3 * t4 + 0.48509427950676e6 * t6 - 0.18880539857133e6 * t6 * R + 0.24591392413587e5 * t6 * t2 - 0.10716456122957e4 * t6 * t4 + 0.26391807499124e6 * t14 - 0.10270490415437e6 * t14 * R + 0.13380372752977e5 * t14 * t2 - 0.58336891438899e3 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = 0.86053887492359e5 * t23 - 0.33569136205528e5 * t23 * R + 0.43834703412849e4 * t23 * t2 - 0.19152871242158e3 * t23 * t4 + 0.20683287070595e5 * t31 - 0.80997694677585e4 * t31 * R + 0.10613602925936e4 * t31 * t2 - 0.46519515434790e2 * t31 * t4 + 0.39993548477133e4 * t39 - 0.15722478899418e4 * t39 * R + 0.20669884747525e3 * t39 * t2 - 0.90849012237643e1 * t39 * t4;
	t48 = t22 + t47;
	}
	else if (R<7.500000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = 0.11936360247650e6 - 0.45585329671514e5 * R + 0.58066568893482e4 * t2 - 0.24690516822800e3 * t4 + 0.33348662774386e6 * t6 - 0.12607119784185e6 * t6 * R + 0.15938399209522e5 * t6 * t2 - 0.67380684429265e3 * t6 * t4 + 0.18013024726909e6 * t14 - 0.68034078890027e5 * t14 * R + 0.85981899578954e4 * t14 * t2 - 0.36349844105192e3 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = 0.57749436450762e5 * t23 - 0.21856949567627e5 * t23 * R + 0.27679963222639e4 * t23 * t2 - 0.11725404488038e3 * t23 * t4 + 0.13584912161109e5 * t31 - 0.51625108845234e4 * t31 * R + 0.65622117766464e3 * t31 * t2 - 0.27892429690930e2 * t31 * t4 + 0.26020166448580e4 * t39 - 0.99403897841554e3 * t39 * R + 0.12694589416128e3 * t39 * t2 - 0.54180987725473e1 * t39 * t4;
	t48 = t22 + t47;
	}
	else if (R<7.750000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = 0.78440330703712e5 - 0.29216020962400e5 * R + 0.36240823947997e4 * t2 - 0.14990185735918e3 * t4 + 0.22579901600936e6 * t6 - 0.82996153148053e5 * t6 * R + 0.10195059917016e5 * t6 * t2 - 0.41854732018127e3 * t6 * t4 + 0.12121101479350e6 * t14 - 0.44466385899792e5 * t14 * R + 0.54558308925305e4 * t14 * t2 - 0.22383803814682e3 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = 0.38102100480840e5 * t23 - 0.13998015179658e5 * t23 * R + 0.17201384038680e4 * t23 * t2 - 0.70682581840568e2 * t23 * t4 + 0.87270609402482e4 * t31 - 0.32193703961789e4 * t31 * R + 0.39713577921871e3 * t31 * t2 - 0.16377523093333e2 * t31 * t4 + 0.16369171088245e4 * t39 - 0.60799916400212e3 * t39 * R + 0.75473918906161e2 * t39 * t2 - 0.31304554278753e1 * t39 * t4;
	t48 = t22 + t47;
	}
	else if (R<8.000000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = 0.50193146726350e5 - 0.18281627164711e5 * R + 0.22131928725172e4 * t2 - 0.89218437045955e2 * t4 + 0.15055258335394e6 * t6 - 0.53868501797564e5 * t6 * R + 0.64366532911464e4 * t6 * t2 - 0.25689542229440e3 * t6 * t4 + 0.80568477907923e5 * t14 - 0.28733790976343e5 * t14 * R + 0.34258186443438e4 * t14 * t2 - 0.13652568338610e3 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = 0.24884808548952e5 * t23 - 0.88816441092496e4 * t23 * R + 0.10599614915573e4 * t23 * t2 - 0.42287875934729e2 * t23 * t4 + 0.54704318817679e4 * t31 - 0.19587397928962e4 * t31 * R + 0.23447376589190e3 * t31 * t2 - 0.93813074663730e1 * t31 * t4 + 0.95791859120856e3 * t39 - 0.34516102815081e3 * t39 * R + 0.41559320731795e2 * t39 * t2 - 0.16717630332790e1 * t39 * t4;
	t48 = t22 + t47;
	}
	else if (R<8.250000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = 0.31030870220996e5 - 0.11095773475203e5 * R + 0.13149611613288e4 * t2 - 0.51792115746437e2 * t4 + 0.98728199352976e5 * t6 - 0.34434357797203e5 * t6 * R + 0.40073852911012e4 * t6 * t2 - 0.15567592229252e3 * t6 * t4 + 0.52921548875121e5 * t14 - 0.18366192589043e5 * t14 * R + 0.21298688459311e4 * t14 * t2 - 0.82527775118909e2 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = 0.16039266740141e5 * t23 - 0.55645659309456e4 * t23 * R + 0.64532671926928e3 * t23 * t2 - 0.25011427089396e2 * t23 * t4 + 0.34024015514909e4 * t31 - 0.11832284190423e4 * t31 * R + 0.13753484416017e3 * t31 * t2 - 0.53421857275511e1 * t31 * t4 + 0.54891505470453e3 * t39 - 0.19178470196180e3 * t39 * R + 0.22387279958169e2 * t39 * t2 - 0.87292800104453e0 * t39 * t4;
	t48 = t22 + t47;
	}
	else if (R<8.500000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = 0.18263676399253e5 - 0.64531575400235e4 * R + 0.75221983585245e3 * t2 - 0.29055092494869e2 * t4 + 0.63515322875535e5 * t6 - 0.21629675441770e5 * t6 * R + 0.24553025813518e4 * t6 * t2 - 0.92965509777394e2 * t6 * t4 + 0.34364076304732e5 * t14 - 0.11618020745266e5 * t14 * R + 0.13119086224430e4 * t14 * t2 - 0.49478877200194e2 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = 0.10343377231256e5 * t23 - 0.34933333822601e4 * t23 * R + 0.39426822851954e3 * t23 * t2 - 0.14867649685366e2 * t23 * t4 + 0.21543567879581e4 * t31 - 0.72939395957586e3 * t31 * R + 0.82524606649083e2 * t31 * t2 - 0.31195498685173e1 * t31 * t4 + 0.33997671140635e3 * t39 - 0.11580712258063e3 * t39 * R + 0.13177876396815e2 * t39 * t2 - 0.50083088745448e0 * t39 * t4;
	t48 = t22 + t47;
	}
	else if (R<8.750000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = 0.99330864827699e4 - 0.35129493342060e4 * R + 0.40631298810924e3 * t2 - 0.15490118073566e2 * t4 + 0.39858496114745e5 * t6 - 0.13280207173255e5 * t6 * R + 0.14730121968206e4 * t6 * t2 - 0.54444318227154e2 * t6 * t4 + 0.21947245058738e5 * t14 - 0.72356097172675e4 * t14 * R + 0.79633085444321e3 * t14 * t2 - 0.29260141200203e2 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = 0.65421885291196e4 * t23 - 0.21517373697414e4 * t23 * R + 0.23643340351733e3 * t23 * t2 - 0.86780487048872e1 * t23 * t4 + 0.12854457668379e4 * t31 - 0.42271948153345e3 * t31 * R + 0.46445256291151e2 * t31 * t2 - 0.17046733838926e1 * t31 * t4 + 0.20938013791313e3 * t39 - 0.69714214288913e2 * t39 * R + 0.77551813036718e1 * t39 * t2 - 0.28817617791947e0 * t39 * t4;
	t48 = t22 + t47;
	}
	else if (R<9.000000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = 0.46092122706189e4 - 0.16876210328971e4 * R + 0.19770403938823e3 * t2 - 0.75431105032415e1 * t4 + 0.24233852291723e5 * t6 - 0.79231864339338e4 * t6 * R + 0.86078125518389e3 * t6 * t2 - 0.31121234736231e2 * t6 * t4 + 0.13794192019669e5 * t14 - 0.44402772467295e4 * t14 * R + 0.47686428638173e3 * t14 * t2 - 0.17089986226432e2 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = 0.40943752756130e4 * t23 - 0.13124871113963e4 * t23 * R + 0.14051908827788e3 * t23 * t2 - 0.50241700290988e1 * t23 * t4 + 0.84161141709936e3 * t31 - 0.27054770448022e3 * t31 * R + 0.29054196056496e2 * t31 * t2 - 0.10421568035248e1 * t31 * t4 + 0.62422686801033e2 * t39 - 0.19328802479051e2 * t39 * R + 0.19968485254019e1 * t39 * t2 - 0.68811119699659e-1 * t39 * t4;
	t48 = t22 + t47;
	}
	else if (R<9.250000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = 0.13427671872307e4 - 0.59880600510105e3 * R + 0.76724591855332e2 * t2 - 0.30623902242453e1 * t4 + 0.14085348698299e5 * t6 - 0.45403519027923e4 * t6 * R + 0.48491075172373e3 * t6 * t2 - 0.17200104978447e2 * t6 * t4 + 0.84373865713169e4 * t14 - 0.26546754306122e4 * t14 * R + 0.27846408459090e3 * t14 * t2 - 0.97418306045498e1 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = 0.25608465270024e4 * t23 - 0.80131086185939e3 * t23 * R + 0.83721727218225e2 * t23 * t2 - 0.29205640639263e1 * t23 * t4 + 0.45645561092051e3 * t31 - 0.14216243575393e3 * t31 * R + 0.14789166198019e2 * t31 * t2 - 0.51382236432194e0 * t31 * t4 + 0.93978751571521e2 * t39 - 0.29847490735880e2 * t39 * R + 0.31655916650496e1 * t39 * t2 - 0.11209790264957e0 * t39 * t4;
	t48 = t22 + t47;
	}
	else if (R<9.500000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = -0.59430843654249e3 + 0.29434737744314e2 * R + 0.88066737098900e1 * t2 - 0.61489767846355e0 * t4 + 0.76476100494436e4 * t6 - 0.24524366653256e4 * t6 * R + 0.25919018551112e3 * t6 * t2 - 0.90660305203352e1 * t6 * t4 + 0.51203746334984e4 * t14 - 0.15788877751035e4 * t14 * R + 0.16216271642779e3 * t14 * t2 - 0.55507903103837e1 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = 0.15337643968450e4 * t23 - 0.46820314397052e3 * t23 * R + 0.47710082041051e2 * t23 * t2 - 0.16228471206049e1 * t23 * t4 + 0.30152085135432e3 * t31 - 0.91913324543274e2 * t31 * R + 0.93568298509215e1 * t31 * t2 - 0.31806249595805e0 * t31 * t4 + 0.33706907326330e2 * t39 - 0.10299865575277e2 * t39 * R + 0.10523348909304e1 * t39 * t2 - 0.35944505384015e-1 * t39 * t4;
	t48 = t22 + t47;
	}
	else if (R<9.750000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = -0.16097376230044e4 + 0.35009658610069e3 * R - 0.24947205064465e2 * t2 + 0.56944894519803e0 * t4 + 0.34786467420237e4 * t6 - 0.11359219366667e4 * t6 * R + 0.12060968775755e3 * t6 * t2 - 0.42035569149468e1 * t6 * t4 + 0.28047275169741e4 * t14 - 0.84763079093796e3 * t14 * R + 0.85188297041947e2 * t14 * t2 - 0.28499334898278e1 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = 0.88181974248145e3 * t23 - 0.26232588469781e3 * t23 * R + 0.26038791591293e2 * t23 * t2 - 0.86245096447301e0 * t23 * t4 + 0.18077605919321e3 * t31 - 0.53783390176607e2 * t31 * R + 0.53431525491672e1 * t31 * t2 - 0.17723171344036e0 * t31 * t4 - 0.11985796425607e2 * t39 + 0.41294092937553e1 * t39 * R - 0.46653614791513e0 * t39 * t2 + 0.17349215277230e-1 * t39 * t4;
	t48 = t22 + t47;
	}
	else if (R<10.000000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = -0.22494763997206e4 + 0.54693928662878e3 * R - 0.45136199990423e2 * t2 + 0.12596709939487e1 * t4 + 0.18585173994554e4 * t6 - 0.63742060049186e3 * t6 * R + 0.69481345585769e2 * t6 * t2 - 0.24555794048003e1 * t6 * t4 + 0.20449555794185e4 * t14 - 0.61385481015160e3 * t14 * R + 0.61211273371552e2 * t14 * t2 - 0.20302061848570e1 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = 0.68765599233129e3 * t23 - 0.20258319234393e3 * t23 * R + 0.19911335965253e2 * t23 * t2 - 0.65296530204428e0 * t23 * t4 + 0.11158271047886e3 * t31 - 0.32493129033730e2 * t31 * R + 0.31595360216928e1 * t31 * t2 - 0.10257815694551e0 * t31 * t4 + 0.46101752969382e2 * t39 - 0.13743682827779e2 * t39 * R + 0.13666015055756e1 * t39 * t2 - 0.45322157491684e-1 * t39 * t4;
	t48 = t22 + t47;
	}
	else if (R<11.000000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = -0.21921515764673e4 + 0.52974183965278e3 * R - 0.43416455292823e2 * t2 + 0.12023461706954e1 * t4 - 0.94668590595597e3 * t6 + 0.20414039113156e3 * t6 * R - 0.14674753576572e2 * t6 * t2 + 0.34962390061099e0 * t6 * t4 + 0.28109886950547e3 * t14 - 0.84697797177699e2 * t14 * R + 0.82955720741618e1 * t14 * t2 - 0.26634947494402e0 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = 0.15181224305327e3 * t23 - 0.41830067560524e2 * t23 * R + 0.38360234869121e1 * t23 * t2 - 0.11712155276625e0 * t23 * t4 + 0.38196342550554e2 * t31 - 0.10477218655238e2 * t31 * R + 0.95794498384341e0 * t31 * t2 - 0.29191789017201e-1 * t31 * t4 + 0.26352916384051e1 * t39 - 0.70374442848669e0 * t39 * R + 0.62607665646337e-1 * t39 * t2 - 0.18556961607085e-2 * t39 * t4;
	t48 = t22 + t47;
	}
	else if (R<12.000000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = -0.15203208371270e4 + 0.34651527437816e3 * R - 0.26759494813311e2 * t2 + 0.69758979252839e0 * t4 - 0.13208814822433e4 * t6 + 0.30619373011899e3 * t6 * R - 0.23952329848158e2 * t6 * t2 + 0.63076257550751e0 * t6 * t4 - 0.22018867312808e3 * t14 + 0.52016987176906e2 * t14 * R - 0.41330446853477e1 * t14 * t2 + 0.11027527534415e0 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = -0.15300701795912e2 * t23 + 0.37461901256181e1 * t23 * R - 0.30727266637335e0 * t23 * t2 + 0.84328761211879e-2 * t23 * t4 - 0.36833808736200e1 * t31 + 0.94452409680975e0 * t31 * R - 0.80395266342738e-1 * t31 * t2 + 0.22730670490456e-2 * t31 * t4 + 0.12164586240584e1 * t39 - 0.31678997002851e0 * t39 * R + 0.27429987604686e-1 * t39 * t2 - 0.78970591702199e-3 * t39 * t4;
	t48 = t22 + t47;
	}
	else if (R<14.000000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = -0.68431251874209e3 + 0.13751319478192e3 * R - 0.93426548469589e1 * t2 + 0.21378868235193e0 * t4 - 0.52500737914619e3 * t6 + 0.10722520434473e3 * t6 * R - 0.73716193669688e1 * t6 * t2 + 0.17018728436338e0 * t6 * t4 - 0.69249641481517e2 * t14 + 0.14282229265264e2 * t14 * R - 0.98848152604417e0 * t14 * t2 + 0.22926298696833e-1 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = -0.13896955168161e1 * t23 + 0.26843855584410e0 * t23 * R - 0.17460035558854e-1 * t23 * t2 + 0.38252526522966e-3 * t23 * t4 + 0.56048874000816e0 * t31 - 0.11644330659725e0 * t31 * R + 0.80186839411803e-2 * t31 * t2 - 0.18287601439653e-3 * t31 * t4 - 0.68567656470026e0 * t39 + 0.15874382716116e0 * t39 * R - 0.12197828827787e-1 * t39 * t2 + 0.31106676165780e-3 * t39 * t4;
	t48 = t22 + t47;
	}
	else if (R<16.000000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = -0.20056145805663e3 + 0.33852253206467e2 * R - 0.19383018772835e1 * t2 + 0.37494564026326e-1 * t4 - 0.12576709010588e3 * t6 + 0.21673713836092e2 * t6 * R - 0.12607986163519e1 * t6 * t2 + 0.24691552205833e-1 * t6 * t4 - 0.15120945142735e2 * t14 + 0.26832229069534e1 * t14 * R - 0.15998107187912e0 * t14 * t2 + 0.32000974071901e-2 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = -0.68829535087447e0 * t23 + 0.11813852028517e0 * t23 * R - 0.67243187332162e-2 * t23 * t2 + 0.12691295985733e-3 * t23 * t4 + 0.46139981187220e0 * t31 - 0.95209964853830e-1 * t31 * R + 0.65020166737928e-2 * t31 * t2 - 0.14676488898255e-3 * t31 * t4 + 0.93328996866654e0 * t39 - 0.18817757284602e0 * t39 * R + 0.12582271172726e-1 * t39 * t2 - 0.27893561930678e-3 * t39 * t4;
	t48 = t22 + t47;
	}
	else if (R<18.000000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = -0.10983164798790e3 + 0.16840413818581e2 * R - 0.87506191554061e0 * t2 + 0.15343731490016e-1 * t4 - 0.62497668472459e2 * t6 + 0.98106972798243e1 * t6 * R - 0.51936008158522e0 * t6 * t2 + 0.92449160648604e-2 * t6 * t4 - 0.54096700430281e1 * t14 + 0.86235882575848e0 * t14 * R - 0.46177066804432e-1 * t14 * t2 + 0.82918063480075e-3 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = -0.81995590162102e0 * t23 + 0.14282487355015e0 * t23 * R - 0.82672158122775e-2 * t23 * t2 + 0.15905664900444e-3 * t23 * t4 - 0.63623951473627e0 * t31 + 0.11059740888526e0 * t31 * R - 0.63609441849000e-2 * t31 * t2 + 0.12121346224022e-3 * t31 * t4 - 0.10406221700827e1 * t39 + 0.18193095316946e0 * t39 * R - 0.10549511703242e-1 * t39 * t2 + 0.20297652394253e-3 * t39 * t4;
	t48 = t22 + t47;
	}
	else if (R<20.000000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = -0.47483106986854e2 + 0.64489903184057e1 * R - 0.29776060997530e0 * t2 + 0.46529665721400e-2 * t4 - 0.20985769209011e2 * t6 + 0.28920474025831e1 * t6 * R - 0.13499064396071e0 * t6 * t2 + 0.21269635162582e-2 * t6 * t4 - 0.15240452466147e1 * t14 + 0.21475469302291e0 * t14 * R - 0.10199059430234e-1 * t14 * t2 + 0.16292123898225e-3 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = 0.52606600127203e0 * t23 - 0.81512110265363e-1 * t23 * R + 0.41959499552507e-2 * t23 * t2 - 0.71742717060895e-4 * t23 * t4 + 0.20286370401257e0 * t31 - 0.29253127572888e-1 * t31 * R + 0.14085300627746e-2 * t31 * t2 - 0.22665690494494e-4 * t31 * t4 + 0.58371166779617e0 * t39 - 0.88791353143683e-1 * t39 * R + 0.44906164252660e-2 * t39 * t2 - 0.75544367326128e-4 * t39 * t4;
	t48 = t22 + t47;
	}
	else if (R<25.000000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = -0.19253689207769e2 + 0.22145776515429e1 * R - 0.86039976632157e-1 * t2 + 0.11242893497543e-2 * t4 - 0.76874142674691e1 * t6 + 0.89729416135182e0 * t6 * R - 0.35252981899145e-1 * t6 * t2 + 0.46466914856546e-3 * t6 * t4 - 0.44563515976596e0 * t14 + 0.52993179995595e-1 * t14 * R - 0.21109837788679e-2 * t14 * t2 + 0.28119978126156e-4 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = -0.12821495282085e0 * t23 + 0.16630032848568e-1 * t23 * R - 0.71115720044580e-3 * t23 * t2 + 0.10042402200714e-4 * t23 * t4 + 0.56801685577257e-1 * t31 - 0.73438248075890e-2 * t31 * R + 0.31306492450971e-3 * t31 * t2 - 0.44079381900798e-5 * t31 * t4 - 0.45076727992929e-1 * t39 + 0.55269062246847e-2 * t39 * R - 0.22529654315223e-3 * t39 * t2 + 0.30541821475108e-5 * t39 * t4;
	t48 = t22 + t47;
	}
	else if (R<30.000000)
	{
	t2 = R * R;
	t4 = t2 * R;
	t6 = legP[2];
	t14 = legP[4];
	t22 = -0.23268011871778e1 + 0.18335108907194e0 * R - 0.47909141333223e-2 * t2 + 0.40968516436477e-4 * t4 - 0.42379419711930e0 * t6 + 0.25659752909839e-1 * t6 * R - 0.38760556146587e-3 * t6 * t2 - 0.20253593691940e-6 * t6 * t4 + 0.13695588670909e-1 * t14 - 0.21265098168275e-2 * t14 * R + 0.93803813629015e-4 * t14 * t2 - 0.12771897738033e-5 * t14 * t4;
	t23 = legP[6];
	t31 = legP[8];
	t39 = legP[10];
	t47 = 0.95414732243320e-1 * t23 - 0.10205529359133e-1 * t23 * R + 0.36226528786222e-3 * t23 * t2 - 0.42698976433936e-5 * t23 * t4 - 0.39348855957060e-1 * t31 + 0.41942401765291e-2 * t31 * R - 0.14845767485500e-3 * t31 * t2 + 0.17456964681165e-5 * t31 * t4 + 0.82329864310103e-2 * t39 - 0.87025950618816e-3 * t39 * R + 0.30590086082683e-4 * t39 * t2 - 0.35763957562138e-6 * t39 * t4;
	t48 = t22 + t47;
	}
	else
	{
	t2 = legP[2];
	t6 = legP[4];
	t10 = legP[6];
	t14 = legP[8];
	t18 = legP[10];
	t22 = 0.89650383248660e16 * R - 0.25974712800427e19 + 0.14110167093330e17 * t2 * R - 0.10294747413919e19 * t2 + 0.15474570148836e16 * t6 * R - 0.58116685316478e17 * t6 + 0.12678113196716e15 * t10 * R - 0.37761566962528e16 * t10 + 0.15270423711763e14 * t14 * R - 0.44013495545803e15 * t14 - 0.34930749728818e14 * t18 * R + 0.10486644683026e16 * t18;
	t23 = R * R;
	t24 = t23 * t23;
	t48 = 0.10000000000000e-10 * t22 / t24 / t23;
	}

	free(legP);

	return t48;
}
