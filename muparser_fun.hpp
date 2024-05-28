#include <muParser.h>
#include <cmath>
#include <memory>
#include <string>

class MuparserFun
{
public:
    // Copy constructor
    MuparserFun(const MuparserFun &m)
        : m_var(0), m_var2(0), m_parser(m.m_parser)
    {
        try
        {
            m_parser.DefineVar("x", &m_var);
            m_parser.DefineVar("y", &m_var2);
            m_parser.DefineConst("p", M_PI);
        }
        catch (mu::Parser::exception_type &e)
        {
            std::cerr << "Error in copy constructor: " << e.GetMsg() << std::endl;
        }
    };

    // Constructor taking an expression string
    MuparserFun(const std::string &s)
        : m_var(0), m_var2(0)
    {
        try
        {
            std::cout << "Defining variables x and y" << std::endl;
            m_parser.DefineVar("x", &m_var);
            m_parser.DefineVar("y", &m_var2);

            std::cout << "Defining constant pi" << std::endl;
            m_parser.DefineConst("p", M_PI); // Use M_PI from cmath

            std::cout << "Setting expression: " << s << std::endl;
            m_parser.SetExpr(s);
        }
        catch (mu::Parser::exception_type &e)
        {
            std::cerr << "Error in constructor: " << e.GetMsg() << std::endl;
        }
    };

    // Function call operator to evaluate the expression
    double operator()(const double &x, const double &y)
    {
        m_var = x;
        m_var2 = y;
        try
        {
            return m_parser.Eval();
        }
        catch (mu::Parser::exception_type &e)
        {
            std::cerr << "Error during evaluation: " << e.GetMsg() << std::endl;
            return std::numeric_limits<double>::quiet_NaN(); // Return NaN in case of error
        }
    };

private:
    double m_var;
    double m_var2;
    mu::Parser m_parser;

};