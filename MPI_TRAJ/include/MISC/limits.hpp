#pragma once
#include "limit.hpp"

#include <string>
#include <vector>

struct Limits
{
	std::vector<Limit> limits; 

	Limits* add_limit( int var_number, double lb, double ub )
	{
		limits.emplace_back( var_number, lb, ub );
		return this;
	}

	Limits* add_limit( int var_number, double lb, std::string ub_str )
	{
		limits.emplace_back( var_number, lb, ub_str );
		return this;
	}

	Limits* add_limit( int var_number, std::string lb_str, double ub )
	{
		limits.emplace_back( var_number, lb_str, ub );
		return this;
	}

	Limits* add_limit( int var_number, std::string lb_str, std::string ub_str )
	{
		limits.emplace_back( var_number, lb_str, ub_str );
		return this;
	}

	void show_limits( void )
	{
		for ( auto& l : limits )
			l.show_limit();
	}

	Limits( ) { } 
};
