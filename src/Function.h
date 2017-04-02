#ifndef FUNCTION_INCLUDED
#define FUNCTION_INCLUDED

class Function
{
public:
	virtual float eval(const float pos[3])=0;
};

#endif