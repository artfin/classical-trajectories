#include <cmath>
#include "co2_ar_dipole.h"

double dipx(double R, double Theta)
{
	double t1,t2,t3,t4,t5,dip;

	t1 = cos(Theta);
	t2 = sin(Theta);
	t3 = t1*t1;
	t4 = t2*t2;
	t1 = t2 * t1;
	t2 = 1 / R;
	t5 = t2*t2;
	dip = t5 * t5 * (t1 * (0.105086403000e3 + (0.293481514286e2 + (-0.377666170574e4 * t3 - 0.859715365822e3 * t4
	 	- 0.189568894120e4) * t2) * t5) +  t5 * (t2 * (0.793818108945e1 * t1 * (-0.178506600000e3 * t3 - 0.429246000000e2 * t4) 
		- 0.793818108945e1 * t1 * (-0.108839000000e3 * t3
	 	+ 0.267430000000e2 * t4) + 0.158763621789e2 * t1 * (-0.108839000000e3 * t4 + 0.267430000000e2 * t3)) +
	 	 0.189771428571e2 * t1 * (0.463950000000e1 * t3 - 0.618600000000e1 * t4) - 0.253028571429e2 * 
	  	t1 * (0.463950000000e1 * t4 - 0.618600000000e1 * t3)));

	return dip;
}


double dipz(double R, double Theta)
{
	double t1, t2, t3,t4,t5, t6, t7, t8, t9, t10, t11, dip;
	t1 = cos(Theta);
	t2 = sin(Theta);

	t2 = t2*t2;
	t1 = t1*t1;
	t3 = t1*t1;
	t4 = 0.158215000000e1 * t1 - 0.31643e1 * t2;
	t5 = 0.158215000000e1 * t2 - 0.31643e1 * t1;
	t6 = 0.128073000000e2 * t2 + 0.266861000000e2 * t1;
	t7 = 0.742320000000e1 * t1;
	t8 = t2*t2 + t3;
	t9 = 1 / R;
	t10 = t9*t9;
	t11 = t9 * t10;
	dip = t10*t10 * ((0.442800000000e2 * t5 - 0.221400000000e2 * t4 - 0.350288010000e2) * t11 * t6 + 
		0.793818108945e1 * (-0.168811200000e3 * t2 - 0.346374700000e3 * t1 - 0.199239900000e3 * t8 -
		0.262897800000e3 * t2 * t1 - 0.147134800000e3 * t3 - (0.360184000000e3 * t1 + 0.774672000000e2 * t2) * t2) *
		t11 + t10 * (0.234785211428e2 * t2 - 0.498918574286e2 * t1 - 0.469570422857e2 * t8 + 0.316960035429e3 * t2 *
		t1 + (-0.145847316995e4 * t2 * t1 - 0.467193074784e4 * t3 + 0.158763621789e2 * (-0.572357600000e3 * t1 -
		0.142506000000e3 * t2) * t2 + 0.287667635578e3 * t6) * t9 - 0.357069175714e2 * t3 + 0.474428571429e1 *
		(t7 - 0.247440000000e1 * t2) * t2 - 0.440222271429e1 + 0.126514285714e2 * (t7 - 0.927900000000e0 * t2) * t2) + 
		0.221400000000e2 * t5 - 0.110700000000e2 * t4 - 0.175144005000e2);

	return dip;
}



double dipx_cos(double R_inv,double cos_theta)
{
	double t1,t2,t3,t4,t5,dip;

	t1 = cos_theta;
	t2 = sqrt(1-t1*t1);
	t3 = t1*t1;
	t4 = 1-t1*t1;
	/*t1 = cos(cos_theta);
	t2 = sin(cos_theta);
	t3 = t1*t1;
	t4 = t2*t2;*/


	t1 = t2 * t1;
	t2 = 1/R_inv;
	t5 = t2*t2;
	dip = t5 * t5 * (t1 * (0.105086403000e3 + (0.293481514286e2 + (-0.377666170574e4 * t3 - 0.859715365822e3 * t4
	 	- 0.189568894120e4) * t2) * t5) +  t5 * (t2 * (0.793818108945e1 * t1 * (-0.178506600000e3 * t3 - 0.429246000000e2 * t4) 
		- 0.793818108945e1 * t1 * (-0.108839000000e3 * t3
	 	+ 0.267430000000e2 * t4) + 0.158763621789e2 * t1 * (-0.108839000000e3 * t4 + 0.267430000000e2 * t3)) +
	 	 0.189771428571e2 * t1 * (0.463950000000e1 * t3 - 0.618600000000e1 * t4) - 0.253028571429e2 * 
	  	t1 * (0.463950000000e1 * t4 - 0.618600000000e1 * t3)));


	return dip;
}

double dipz_cos(double R_inv, double cos_theta)
{
	double t1, t2, t3,t4,t5, t6, t7, t8, t9, t10, t11, dip;
	t1 = cos_theta;
	//t2 = sqrt(1-t1*t1);
	t2 = 1-t1*t1;
	t1 = t1*t1;
	t3 = t1*t1;

	/*t1 = cos(cos_theta);
	t2 = sin(cos_theta);
	t2 = t2*t2;
	t1 = t1*t1;
	t3 = t1*t1;*/


	t4 = 0.158215000000e1 * t1 - 0.31643e1 * t2;
	t5 = 0.158215000000e1 * t2 - 0.31643e1 * t1;
	t6 = 0.128073000000e2 * t2 + 0.266861000000e2 * t1;
	t7 = 0.742320000000e1 * t1;
	t8 = t2*t2 + t3;
	t9 = 1/R_inv;
	t10 = t9*t9;
	t11 = t9 * t10;
	dip = t10*t10 * ((0.442800000000e2 * t5 - 0.221400000000e2 * t4 - 0.350288010000e2) * t11 * t6 + 
		0.793818108945e1 * (-0.168811200000e3 * t2 - 0.346374700000e3 * t1 - 0.199239900000e3 * t8 -
		0.262897800000e3 * t2 * t1 - 0.147134800000e3 * t3 - (0.360184000000e3 * t1 + 0.774672000000e2 * t2) * t2) *
		t11 + t10 * (0.234785211428e2 * t2 - 0.498918574286e2 * t1 - 0.469570422857e2 * t8 + 0.316960035429e3 * t2 *
		t1 + (-0.145847316995e4 * t2 * t1 - 0.467193074784e4 * t3 + 0.158763621789e2 * (-0.572357600000e3 * t1 -
		0.142506000000e3 * t2) * t2 + 0.287667635578e3 * t6) * t9 - 0.357069175714e2 * t3 + 0.474428571429e1 *
		(t7 - 0.247440000000e1 * t2) * t2 - 0.440222271429e1 + 0.126514285714e2 * (t7 - 0.927900000000e0 * t2) * t2) + 
		0.221400000000e2 * t5 - 0.110700000000e2 * t4 - 0.175144005000e2);


	return dip;
}


double ddipxdR(double R, double Theta)
{
	double t1, t2, t3, t4, t5, t10, t12, t16, t17, t37, t41, t44;

	t1 = cos(Theta);
	t2 = sin(Theta);
	t3 = t1 * t2;
	t4 = R * R;
	t5 = t4 * t4;
	t10 = t1 * t1;
	t12 = t2 * t2;
	t16 = t5 * t5;
	t17 = 0.1e1 / t16;
	t37 = 0.1e1 / t5 / t4 / R;
	t41 = t10 * t1 * t2;
	t44 = t12 * t2 * t1;
    return  -0.420345612000e3 * t3 / t5 / R + 0.735604821000e3 * (0.128073000000e2 * t10 + 0.266861000000e2 * t12) * t1 * t2 * t17 +
    		 0.215093642400e4 * t3 * (0.158215000000e1 * t10 - 0.316430000000e1 * t12) * t17 + 0.132698225884e5 * t3 * t17 - 0.430187284800e4 * 
    		 t3 * (0.158215000000e1 * t12 - 0.316430000000e1 * t10) * t17 - 0.176088908572e3 * t3 * t37 - 0.113862857143e3 * (0.463950000000e1 * 
    		 t41 - 0.618600000000e1 * t44) * t37 + 0.151817142857e3 * (0.463950000000e1 * t44 - 0.618600000000e1 * t41) * t37 - 0.555672676262e2 *
    		 (-0.178506600000e3 * t41 - 0.429246000000e2 * t44) * t17 + 0.555672676262e2 * (-0.108839000000e3 * t41 + 0.267430000000e2 * t44) * 
    		 t17 - 0.111134535252e3 * (-0.108839000000e3 * t44 + 0.267430000000e2 * t41) * t17;

}

double ddipxdTheta(double R, double Theta)
{
	double t1, t2 , t3, t4, t5, t8, t9, t12, t15, t20, t29, t42, t50, t55, t56, t58, t68;

	t1 = sin(Theta);
	t2 = t1 * t1;
	t3 = R * R;
	t4 = t3 * t3;
	t5 = 0.1e1 / t4;
	t8 = cos(Theta);
	t9 = t8 * t8;
	t12 = t9 * t2;
	t15 = 0.1e1 / t4 / t3 / R;
	t20 = 0.128073000000e2 * t9 + 0.266861000000e2 * t2;
	t29 = 0.158215000000e1 * t9 - 0.316430000000e1 * t2;
	t42 = 0.158215000000e1 * t2 - 0.316430000000e1 * t9;
	t50 = 0.1e1 / t4 / t3;
	t55 = 0.324765000000e2 * t12;
	t56 = t9 * t9;
	t58 = t2 * t2;
	t68 = 0.406746000000e3 * t12;
	double ddip = -0.105086403000e3 * t2 * t5 + 0.105086403000e3 * t9 * t5 + 0.583389267983e4 * t12 * t15 + 0.105086403000e3 * t20 * t2 * 
				t15 - 0.105086403000e3 * t20 * t9 * t15 + 0.307276632000e3 * t2 * t29 * t15 - 0.307276632000e3 * t9 * t29 * t15 + 0.189568894120e4 * 
				t2 * t15 - 0.189568894120e4 * t9 * t15 - 0.614553264000e3 * t2 * t42 * t15 + 0.614553264000e3 * t9 * t42 * t15 - 0.293481514286e2 *
				 t2 * t50 + 0.293481514286e2 * t9 * t50 + 0.189771428571e2 * (-t55 + 0.463950000000e1 * t56 + 0.618600000000e1 * t58) * t50 -
				  0.253028571429e2 * (t55 - 0.463950000000e1 * t58 - 0.618600000000e1 * t56) * t50 + 0.793818108945e1 * (t68 - 0.178506600000e3 * 
				  	t56 + 0.429246000000e2 * t58) * t15 - 0.793818108945e1 * (t68 - 0.108839000000e3 * t56 - 0.267430000000e2 * t58) * t15 + 
				  0.158763621789e2 * (-t68 + 0.108839000000e3 * t58 + 0.267430000000e2 * t56) * t15;
	return ddip;

}


double ddipzdR(double R, double Theta)
{
	double t1, t2, t3, t4, t5, t8, t9, t11, t15, t20, t21, t22, t27, t40, t42, t43, t48;

	t1 = cos(Theta);
	t2 = t1 * t1;
	t3 = 0.158215000000e1 * t2;
	t4 = sin(Theta);
	t5 = t4 * t4;
	t8 = R * R;
	t9 = t8 * t8;
	t11 = 0.1e1 / t9 / R;
	t15 = 0.158215000000e1 * t5;
	t20 = t2 * t5;
	t21 = t9 * t9;
	t22 = 0.1e1 / t21;
	t27 = 0.128073000000e2 * t5 + 0.266861000000e2 * t2;
	t40 = t2 * t2;
	t42 = 0.742320000000e1 * t20;
	t43 = t5 * t5;
	t48 = 0.1e1 / t9 / t8 / R;
	double ddip = 0.442800000000e2 * (t3 - 0.31643e1 * t5) * t11 + 0.700576020000e2 * t11 - 0.885600000000e2 * (t15 - 0.31643e1 * t2) * 
		 t11 + 0.102093121896e5 * t20 * t22 + 0.154980000000e3 * t27 * (t3 - 0.316430000000e1 * t5) * t22 + 0.245201607000e3 * t27 * 
		 t22 - 0.309960000000e3 * t27 * (t15 - 0.316430000000e1 * t2) * t22 - 0.284657142857e2 * (-0.927900000000e0 * t40 + t42 - 
		 0.247440000000e1 * t43) * t48 - 0.569314285714e2 * (-0.309300000000e0 * t2 + 0.123720000000e1 * t5) * t48 + 0.227725714286e3 * 
		 (-0.835110000000e1 * t20 + 0.123720000000e1 * t43 + 0.123720000000e1 * t40) * t48 + 0.264133362857e2 * t48 + 0.227725714286e3 * 
		 (-0.309300000000e0 * t5 + 0.123720000000e1 * t2) * t48 - 0.759085714284e2 * (-0.927900000000e0 * t43 + t42 - 0.247440000000e1 * t40) * 
		 t48 - 0.555672676262e2 * (-0.262897800000e3 * t20 - 0.199239900000e3 * t43 - 0.199239900000e3 * t40) * t22 - 0.555672676262e2 * 
		 (-0.103772400000e3 * t5 - 0.199239900000e3 * t2) * t22 + 0.555672676262e2 * (0.360184000000e3 * t20 + 0.147134800000e3 * t40 + 
		 	0.774672000000e2 * t43) * t22 + 0.555672676262e2 * (0.650388000000e2 * t5 + 0.147134800000e3 * t2) * t22 - 0.111134535252e3 * 
		 (-0.142506000000e3 * t43 - 0.572357600000e3 * t20 - 0.294269600000e3 * t40) * t22 - 0.201367344905e4 * (0.128073e2 * t5 + 0.266861e2 * t2) * t22;

	return ddip;
}


double ddipzdTheta(double R, double Theta)
{
	double t1, t2, t3, t4, t5, t9, t11, t14, t17, t19, t47;

	t1 = cos(Theta);
	t2 = sin(Theta);
	t3 = t1 * t2;
	t4 = R * R;
	t5 = t4 * t4;
	t9 = t2 * t2;
	t11 = t1 * t9 * t2;
	t14 = 0.1e1 / t5 / t4 / R;
	t17 = t1 * t1;
	t19 = t17 * t1 * t2;
	t47 = 0.1e1 / t5 / t4;
	double ddip = 0.315259209000e3 * t3 / t5 + 0.291694633990e4 * t11 * t14 - 0.291694633990e4 * t19 * t14 + 0.614553264000e3 * 
			t3 * (0.158215000000e1 * t17 - 0.316430000000e1 * t9) * t14 + 0.630518418000e3 * (0.128073000000e2 * t9 + 
			0.266861000000e2 * t17) * t1 * t2 * t14 - 0.419358527893e4 * t3 * t14 - 0.122910652800e4 * t3 * (0.158215000000e1 * 
			t9 - 0.316430000000e1 * t17) * t14 + 0.474428571429e1 * (0.185580000000e2 * t19 - 0.247440000000e2 * t11) * t47 + 
			0.146740757143e3 * t3 * t47 - 0.379542857143e2 * (0.216510000000e2 * t11 - 0.216510000000e2 * t19) * t47 + 0.126514285714e2 * 
			(-0.185580000000e2 * t11 + 0.247440000000e2 * t19) * t47 + 0.793818108945e1 * (-0.271164000000e3 * t11 + 0.271164000000e3 * t19) * 
			t14 - 0.793818108945e1 * (-0.410499200000e3 * t11 + 0.131828800000e3 * t19) * t14 + 0.158763621789e2 * (0.574691200000e3 * t11 + 
			0.3236320000e2 * t19) * t14;


	return ddip;
}