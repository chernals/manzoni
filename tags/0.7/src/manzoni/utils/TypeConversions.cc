
template <class T> inline T
s2n(std::string const& s)
{
    std::istringstream i(s);
    T x;
    if(!(i>>x))
    {
        throw("Conversion error");
    }
    return x;
}
  
inline std::string
i2s(int const& i)
{
    std::string s;
    std::stringstream out;
    out << i;
    s = out.str();
    return s;
}

inline std::string
d2s(double const& i)
{
    std::string s;
    std::stringstream out;
    out << i;
    s = out.str();
    return s;
}

