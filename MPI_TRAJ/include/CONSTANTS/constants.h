#pragma once

#include <cmath>

namespace constants
{
	// hartree to joules
	const double HTOJ = 4.35974417 * pow(10, -18);

	// boltzmann constant
	const double BOLTZCONST = 1.38064852 * pow(10, -23);

	// atomic mass unit == m_e ( electron mass )
	const double AMU = 9.10938356 * 1E-31; 

	// atomic length unit to meters == a0
	const double ALU = 0.52917721067 * 1E-10;  

	// atomic time unit to seconds
	const double ATU = 2.418884326505 * pow( 10, -17 );

	// atomic unit of velocity == a0 * Eh / hbar
	const double AVU = 2.18769126277 * 1E6;

	// convert herz to cm^-1
	const double HZTOCM = 3.335631 * pow( 10, -11 );	

	// vacuum permittivity
	const double EPSILON0 = 8.854187817 * 1E-12;
	
	// atomic unit of dipole moment == e * a0
	const double ADIPMOMU = 8.478353552 * 1E-30; 
}
