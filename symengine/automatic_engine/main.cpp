#include <chrono>

#include <symengine/lambda_double.h>
#include <symengine/symengine_exception.h>
//#include <symengine/eval_mpfr.h>
#include <symengine/visitor.h>

#ifdef HAVE_SYMENGINE_LLVM
#include <symengine/llvm_double.h>
using SymEngine::LLVMDoubleVisitor;
#endif

using SymEngine::Basic;
using SymEngine::RCP;
using SymEngine::real_double;
using SymEngine::symbol;
using SymEngine::add;
using SymEngine::mul;
using SymEngine::pow;
using SymEngine::integer;
using SymEngine::vec_basic;
using SymEngine::complex_double;
using SymEngine::LambdaRealDoubleVisitor;
using SymEngine::LambdaComplexDoubleVisitor;
using SymEngine::max;
using SymEngine::sin;
using SymEngine::cos;
using SymEngine::tan;
using SymEngine::cot;
using SymEngine::csc;
using SymEngine::sec;
using SymEngine::asin;
using SymEngine::acos;
using SymEngine::atan;
using SymEngine::acot;
using SymEngine::acsc;
using SymEngine::asec;
using SymEngine::sinh;
using SymEngine::cosh;
using SymEngine::tanh;
using SymEngine::coth;
using SymEngine::csch;
using SymEngine::sech;
using SymEngine::asinh;
using SymEngine::acosh;
using SymEngine::atanh;
using SymEngine::acoth;
using SymEngine::acsch;
using SymEngine::asech;
using SymEngine::atan2;
using SymEngine::log;
using SymEngine::E;
using SymEngine::gamma;
using SymEngine::loggamma;
using SymEngine::min;
using SymEngine::NotImplementedError;
using SymEngine::SymEngineException;

int main(int argc, char *argv[])
{
    //SymEngine::print_stack_on_segfault();

    RCP<const Basic> q1, q2, mu1, mu2, l, r;
    q1 = symbol("q1");
    q2 = symbol("q2");
    mu1 = symbol("mu1");
    mu2 = symbol("mu2");
    l = symbol("l");

    r = add(q1, add(mul(q1, q2), pow(q1, integer(2))));

    LambdaRealDoubleVisitor v;
    v.init({q1, q2}, *r);

    //auto t1 = std::chrono::high_resolution_clock::now();

    //v.call({1.5, 2.0});

    //auto t2 = std::chrono::high_resolution_clock::now();

    //std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << "ms" << std::endl;
    return 0;
}
                












