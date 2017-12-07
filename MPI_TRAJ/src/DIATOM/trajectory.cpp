#include "trajectory.hpp"

Trajectory::Trajectory( Parameters& parameters, const int& N )
{
	this->N = N;
	this->parameters = parameters;
	
	vmblock = vminit();
	y0 = (REAL*) vmalloc(vmblock, VEKTOR, N, 0);
}

void Trajectory::syst (REAL t, REAL *y, REAL *f)
{
  	(void)(t); // avoid unused parameter warning 

	double *out = new double[4];

	rhs( out, y[0], y[1], y[2], y[3] );

	//cout << "out[0]: " << out[0] << endl;
	//cout << "out[1]: " << out[1] << endl;
	//cout << "out[2]: " << out[2] << endl;
	//cout << "out[3]: " << out[3] << endl;

	f[0] = out[0]; // \dot{R} 
	f[1] = out[1]; // \dot{p_R}
	f[2] = out[2]; // \dot{\theta}
	f[3] = out[3]; // \dot{p_\theta}

	delete [] out;
}

//void Trajectory::call_syst( void *self )
//{
	//static_cast< Trajectory* >(self)-kill();
//}

void Trajectory::run_trajectory( void )
{
	// out of memory?
	if ( !vmcomplete(vmblock) )
	{ 
		cout << "mgear: out of memory" << endl;
		return;
	}

	int counter = 0;
	double R_end_value = y0[0] + 0.001;

	t0 = 0.0;
   
  	// initial step size
	h = 0.1;
	xend = parameters.sampling_time; // initial right bound of integration

	while ( y0[0] < R_end_value ) 
	{
		if ( counter == parameters.MaxTrajectoryLength )
		{
			//cout << "Trajectory cut!" << endl;
			break;
		}
		
		//dglsysfnk p = &Trajectory::syst;	
	
		fehler = gear4(&t0, xend, N, std::bind(&Trajectory::syst, this), y0, epsabs, epsrel, &h, fmax, &aufrufe);

		if ( fehler != 0 ) 
		{
			cout << "Gear4: error n = " << 10 + fehler << endl;
			break;
		}
			
		if ( parameters.use_S_matrix == true )
		{
			transform_dipole( temp, y0[0], y0[2] );
		}
		else
		{
			dipole_without_S( temp, y0[0] );
			//cout << "transformed dipole" << endl;
		}
			
		dipx.push_back( temp[0] );
		dipy.push_back( temp[1] );
		dipz.push_back( temp[2] );

		xend = parameters.sampling_time * (counter + 2);

		aufrufe = 0;  // actual number of calls

		counter++;
	}
	
}
