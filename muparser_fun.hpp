/**
 * @file muparser_fun.hpp
 * @author Francesca Visalli
 * @brief Definition of an auxiliary function that uses muParser in order to create a function wrapper from a string.
 * @version 0.1
 * @date 2024-06-02
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#ifndef MUPARSER_FUN_HPP
#define MUPARSER_FUN_HPP

#include <muParser.h>
#include <cmath>
#include <memory>
#include <string>


/**
 * @brief Create a Mu Parser Function object
 * 
 * @param expression A string that contains the expression of the function. Remark:Every p in the string is set as M_PI.
 * @return std::function<double(double, double)> 
 */
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