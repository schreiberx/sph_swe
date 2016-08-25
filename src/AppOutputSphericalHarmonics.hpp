/*
 * OutputSphericalHarmonics.hpp
 *
 *  Created on: 15 Aug 2016
 *      Author: martin
 */

#ifndef SRC_OUTPUTSPHERICALHARMONICS_HPP_
#define SRC_OUTPUTSPHERICALHARMONICS_HPP_


#include <sph/SPHConfig.hpp>
#include <sph/SPHData.hpp>
#include <cassert>

class AppOutputSphericalHarmonics
{
public:
	void run(SPHConfig *sphConfig)
	{
		SPHData h(sphConfig);
/*
		std::cout << std::endl;
		std::cout << "nphi = " << shtns->nphi << std::endl;
		std::cout << "nlat = " << shtns->nlat << std::endl;
		std::cout << std::endl;
		std::cout << "mmax = " << shtns->mmax << std::endl;
		std::cout << "nmax = " << shtns->lmax << std::endl;
		std::cout << std::endl;
		std::cout << "nlm = " << shtns->nlm << std::endl;
*/
		int counter = 0;
		// iterate over modes
		for (int n = 0; n <= sphConfig->spec_n_max; n++)
		{
			for (int m = 0; m <= std::min((int)sphConfig->spec_m_max, n); m++)
			{
				h.spec_update_lambda(
						[&](int i_n, int i_m, cplx &o_data)
						{
							if (i_n == n && i_m == m)
								o_data = 1;
							else
								o_data = 0;
						}
					);

				char buffer[1024];
				sprintf(buffer, "SPH_n%i_m%i.csv", n, m);

				h.spat_write_file(buffer);

				counter++;
			}
		}

		assert(counter == sphConfig->spec_num_elems);
	}


};



#endif /* SRC_OUTPUTSPHERICALHARMONICS_HPP_ */
