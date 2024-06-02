#ifndef MUPARSER_FUN_HPP
#define MUPARSER_FUN_HPP

#include <muParser.h>
#include <cmath>
#include <memory>
#include <string>



inline std::function<double(double, double)> createMuParserFunction(const std::string& expression) {
    return [expression](double x, double y) {
        try {
            mu::Parser p;
            p.DefineConst("p", M_PI);
            p.DefineVar("x", &x);
            p.DefineVar("y", &y);
            p.SetExpr(expression);
            return p.Eval();
        } catch (mu::Parser::exception_type &e) {
            std::cerr << "Error: " << e.GetMsg() << std::endl;
            return 0.0;
        }
    };
}

#endif // MUPARSER_FUN_HPP