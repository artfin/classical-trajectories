#include "spectrum_info.hpp"

void SpectrumInfo::reserve_space( void  )
{
	specfunc_package.reserve( size );
 	zero_out( specfunc_package );

	spectrum_package.reserve( size );
	zero_out( spectrum_package );

	specfunc_chunk.reserve( size );
	zero_out( specfunc_chunk );

	spectrum_chunk.reserve( size );
	zero_out( specfunc_chunk );

	specfunc_total.reserve( size );
	zero_out( specfunc_total );

	spectrum_total.reserve( size );
	zero_out( spectrum_total );
}

void SpectrumInfo::zero_out( std::vector<double>& v )
{
	for ( size_t i = 0; i < size; i++ )
	{
		v[i] = 0;
	}
}

void SpectrumInfo::multiply_chunk( double multiplier )
{
	for ( size_t i = 0; i < size; i++ )
	{
		specfunc_chunk[i] *= multiplier;
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

double SpectrumInfo::receive( int source )
{
	std::cout << "Receiving info in SpectrumInfo object" << std::endl;
	MPI_Status status;
	
	// deriving a source by myself
	if ( source == -1 )
		MPI_Recv( &specfunc_package[0], size, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
	else
	{
		MPI_Recv( &specfunc_package[0], size, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		source = status.MPI_SOURCE;
	}

	MPI_Recv( &spectrum_package[0], size, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
	MPI_Recv( &m2_package, 1, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

	return source;
}
