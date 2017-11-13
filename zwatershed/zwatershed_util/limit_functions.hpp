#pragma once

#include <cstddef>

double LOW_THRESH=  .0001;

struct square
{
private:
    double coef_;

public:
    square( double f )
        : coef_(f)
    {}

    std::size_t operator()( float v ) const
    {
        if ( v < LOW_THRESH ) return 0;
        return v * v * coef_;
    }
    std::size_t operator()(unsigned char v) const
    {
        return (*this)((float)(v) / 255.0f);
    }
};