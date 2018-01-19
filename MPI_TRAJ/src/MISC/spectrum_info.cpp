#include "spectrum_info.hpp"

void SpectrumInfo::reserve_space( void  )
{
	specfunc_package.resize( size );
	spectrum_package.resize( size );
	specfunc_chunk.resize( size );
	spectrum_chunk.resize( size );
	specfunc_total.resize( size );
	spectrum_total.resize( size );
}

void SpectrumInfo::multiply_chunk( double multiplier )
{
	for ( size_t i = 0; i < size; i++ )
	{
		specfunc_chunk[i] *= multiplier;
		//std::cout << "specfunc_chunk[" << i << "]: " << specfunc_chunk[i] << std::endl;

		spectrum_chunk[i] *= multiplier;
	}

	m2_chunk *= multiplier;
}

void SpectrumInfo::add_package_to_chunk( void )
{
	for ( size_t i = 0; i < size; i++ )
	{
		specfunc_chunk[i] += specfunc_package[i];
		spectrum_chunk[i] += spectrum_package[i];
	}

	//for ( int i = 0; i < 5; i++ )
	//{
		//std::cout << "specfunc_chunk[" << i << "] = " << specfunc_chunk[i] << std::endl;
	//}

	m2_chunk += m2_package;
}

void SpectrumInfo::add_chunk_to_total( void )
{
	for ( size_t i = 0; i < size; i++ )
	{
		specfunc_total[i] += specfunc_chunk[i];
		spectrum_total[i] += spectrum_chunk[i];

		specfunc_chunk[i] = 0;
		spectrum_chunk[i] = 0;
	}

	m2_total += m2_chunk;
	m2_chunk = 0;
}

void SpectrumInfo::zero_out_chunk( void )
{
	for ( size_t i = 0; i < size; i++ )
	{
		specfunc_chunk[i] = 0;
		spectrum_chunk[i] = 0;
	}

	m2_chunk = 0;
}

void SpectrumInfo::clear_package( void )
{
	specfunc_package.clear();
	spectrum_package.clear();
	m2_package = 0;
}

void SpectrumInfo::correct( SpectrumInfo& classical, double omega, double freq_step, double kT, bool beware_zero )
{
	double correction;
	
	if ( beware_zero && omega == 0 ) 
		correction = 0.0;
	else 
		correction = corrector( omega, kT );
		
	specfunc_package.push_back( correction * classical.specfunc_package.back() );
	spectrum_package.push_back( correction * classical.spectrum_package.back() );

	m2_package += freq_step * spectrum_package.back();
}

void SpectrumInfo::send( void )
{
	size = specfunc_package.size();
	//std::cout << "(misc) sending specfunc package" << std::endl;
	MPI_Send( &specfunc_package[0], size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
	//std::cout << "(misc) specfunc package sent" << std::endl;
	//std::cout << "(misc) sending spectrum package" << std::endl;
	MPI_Send( &spectrum_package[0], size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
	//std::cout << "(misc) spectrum package" << std::endl;
	//std::cout << "(misc) sending m2 package" << std::endl;
	MPI_Send( &m2_package, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );	
	//std::cout << "(misc) m2 package sent" << std::endl;
}

void SpectrumInfo::receive( int& source, bool define_source )
{
	MPI_Status status;
	
	//std::cout << "(misc) define_source: " << std::boolalpha << define_source << std::endl;
	
	// finding out a source
	if ( define_source )
	{
		MPI_Recv( &specfunc_package[0], size, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		source = status.MPI_SOURCE;
	}
	else
	{
		MPI_Recv( &specfunc_package[0], size, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
	}

	MPI_Recv( &spectrum_package[0], size, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
	MPI_Recv( &m2_package, 1, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
	//std::cout << "(misc) received spectrum package and m2" << std::endl;
}

void SpectrumInfo::save( std::vector<double>& v1, std::vector<double>& v2, std::string filename )
{
	assert( v1.size() == v2.size() && "Sizes of saved vectors should be equal" );

	std::ofstream file( filename );

	for ( size_t k = 0; k < v1.size(); k++ )
	{
		file << v1[k] << " " << v2[k] << endl;
	}
	
	file.close();	
}

void SpectrumInfo::save( const double& m2, std::string filename )
{
	std::ofstream file( filename );

	file << "M2: " << m2 << endl;
		
	file.close();
}

void SpectrumInfo::saving_procedure( Parameters& parameters, std::vector<double>& freqs, std::string filename )
{
	std::string out_dir = parameters.output_directory;
	std::cout << "out_dir: " << out_dir << std::endl;

	bool status = check_dir_exists( out_dir );
	if ( !status )
	{
		mkdir ( out_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
	}

	if ( filename != "" )
	{
		save( freqs, spectrum_chunk, out_dir + "/" + filename );
	}

	std::string specfunc_filename = modify_filename( parameters.specfunc_filename, modifier );
	std::cout << "specfunc_filename: " << specfunc_filename << std::endl;
	save( freqs, specfunc_total, out_dir + "/" + specfunc_filename );

	std::string spectrum_filename = modify_filename( parameters.spectrum_filename, modifier );
	save( freqs, spectrum_total, out_dir + "/" + spectrum_filename );

	std::string m2_filename = modify_filename( parameters.m2_filename, modifier );
	save( m2_total, out_dir + "/" + m2_filename );
}

std::string SpectrumInfo::modify_filename( std::string filename, std::string modifier )
{
	size_t dot_position = filename.find( "." );
	if ( dot_position != std::string::npos )
	{
		std::string str = filename.substr( 0, dot_position );
		str = str + "_" + modifier + ".txt";
		return str;
	}
	else
	{
		std::cout << "There is not dot in filename! Setting filename to base" << std::endl;
		return modifier + ".txt";
	}
}

bool SpectrumInfo::check_dir_exists( std::string dirname )
{
	DIR* dir = opendir( dirname.c_str() );
	
	// directory exists
	if ( dir )
	{
		closedir( dir );
		return true;
	}
	// directory doesn't exist
	else if ( ENOENT == errno ) return false;
	// unknown error
	else
	{
		std::cout << "Unknown error status while checking ouput dir" << std::endl;
		return false;
	}
}



