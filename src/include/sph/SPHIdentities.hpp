/*
 * SPHIdentities.hpp
 *
 *  Created on: 26 Aug 2016
 *      Author: martin
 */

#ifndef SRC_INCLUDE_SPH_SPHIDENTITIES_HPP_
#define SRC_INCLUDE_SPH_SPHIDENTITIES_HPP_

#include <cassert>

class SPHIdentities
{
public:
	inline
	static double D(double k, double m)
	{
		double n=k+1;
		assert(n >= 0);
//		if (n < 0)
//			n = -n-1;
		return ((2.0*n+1.0)*std::sqrt((n*n-m*m)/(4.0*n*n-1.0)));
	}

	inline
	static double E(double n, double m)
	{
		assert(n >= 0);
		return -n;
	}

	inline
	static double R(double k, double m)
	{
		double n=k+1;
//		if (n < 0)
//			n = -n-1;
		if (n < 0)
			return 0;

		assert(n >= 0);
		return std::sqrt((n*n-m*m)/(4.0*n*n-1.0));
	}

	inline
	static double S(double k, double m)
	{
		double n=k-1;
//		if (n < 0)
//			n = -n-1;
		if (n < 0)
			return 0;

		assert(n >= 0);
		return std::sqrt(((n+1.0)*(n+1.0)-m*m)/((2.0*n+1.0)*(2.0*n+3.0)));
	}
};


#endif /* SRC_INCLUDE_SPH_SPHIDENTITIES_HPP_ */
