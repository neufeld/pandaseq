/* PANDAseq -- Assemble paired FASTQ Illumina reads and strip the region between amplification primers.
    Copyright (C) 2011-2012  Andre Masella

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */
#ifndef PROB_H
#        define PROB_H

#        define PROBABILITY(score) (pow(10.0, (-(double)(score)) / 10.0))
#        define PHREDMAX 46
#        define PHREDCLAMP(x) ((x) > PHREDMAX ? PHREDMAX : ((x) < 0 ? 0 : (x)))

#endif
