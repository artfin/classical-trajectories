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
		specfunc_total[i] += specfunc_package[i];
		spectrum_total[i] += spectrum_total[i];
	}

	m2_total += m2_chunk;
}

void SpectrumInfo::add_chunk_to_total_and_zero_out( void )
{
	add_chunk_to_total();
	zero_out_chunk();
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

void SpectrumInfo::zero_out_package( void )
{
	for ( size_t i = 0; i < size; i++ )
	{
		specfunc_package[i] = 0;
		spectrum_package[i] = 0;
	}

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

	MPI_Send( &specfunc_package[0], size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
	MPI_Send( &spectrum_package[0], size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
	MPI_Send( &m2_package, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );	
}

double SpectrumInfo::receive( int source, bool define_source )
{
	//std::cout << "Receiving info in SpectrumInfo object" << std::endl;
	MPI_Status status;
	
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

	return source;
}
