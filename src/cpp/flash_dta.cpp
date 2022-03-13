/* Portions Copyright 2019 Xuesong Zhou and Peiheng Li
 *
 * If you help write or modify the code, please also list your names here.
 * The reason of having Copyright info here is to ensure all the modified version, as a whole, under the GPL
 * and further prevent a violation of the GPL.
 *
 * More about "How to use GNU licenses for your own software"
 * http://www.gnu.org/licenses/gpl-howto.html
 */

#ifdef _WIN32
#include "pch.h"
#endif

#include "config.h"

int main()
{
	generate_default_settings();
	perform_network_assignment();	

	return 0;
}