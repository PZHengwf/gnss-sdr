/*!
 * \file optimize_fft_size.h
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

#ifndef GNSS_SDR_OPTIMIZE_FFT_SIZE_H_
#define GNSS_SDR_OPTIMIZE_FFT_SIZE_H_

#include <iterator>

template <class InputIterator, class T>
InputIterator find_next_larger_or_equal_value(InputIterator first, InputIterator last, const T& val)
{
    while (first != last)
        {
            if (*first >= val) return first;
            ++first;
        }
    return last;
}

unsigned int optimize_fft_size(unsigned int minimum_fft_size);

#endif
