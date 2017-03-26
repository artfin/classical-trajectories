#include <iostream>
#include <chrono>

#include <symengine/basic.h>
#include <symengine/add.h>
#include <symengine/symbol.h>
#include <symengine/dict.h>
#include <symengine/integer.h>
#include <symengine/mul.h>
#include <symengine/pow.h>

using SymEngine::Basic;
using SymEngine::Add;
using SymEngine::Mul;
using SymEngine::Pow;
using SymEngine::Symbol;
using SymEngine::symbol;
using SymEngine::map_vec_int;
using SymEngine::integer;
using SymEngine::multinomial_coefficients;
using SymEngine::RCP;
using SymEngine::rcp_dynamic_cast;

int main(int argc, char *argv[])
{
    SymEngine::print_stack_on_segfault();

    RCP<const Basic> x = symbol("x");
    RCP<const Basic> a, c;

    return 0;
}







