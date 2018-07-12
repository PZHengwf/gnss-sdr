/*!
 * \file optimize_fft_size.cc
 * \brief Computes a FFT size than can be computed fast.
 * \author Carles Fernandez, 2018. cfernandez(at)cttc.es
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 2010-2018  (see AUTHORS file for a list of contributors)
 *
 * GNSS-SDR is a software defined Global Navigation
 *          Satellite Systems receiver
 *
 * This file is part of GNSS-SDR.
 *
 * GNSS-SDR is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GNSS-SDR is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNSS-SDR. If not, see <https://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------------
 */

#include "optimize_fft_size.h"
#include <vector>

unsigned int optimize_fft_size(unsigned int minimum_fft_size)
{
    unsigned int result = 0;
    std::vector<unsigned int> optimum_fft_sizes = {512, 1024, 2000, 2048,
        4096, 6144, 7168, 8192, 10240, 12000, 15000, 16384,
        20000, 30000, 32768, 40000, 50000, 65536, 131072, 262144, 524288,
        1048576, 2097152};
    std::vector<unsigned int>::iterator it;
    it = find_next_larger_or_equal_value(optimum_fft_sizes.begin(),
        optimum_fft_sizes.end(),
        minimum_fft_size);
    return *it;
}
