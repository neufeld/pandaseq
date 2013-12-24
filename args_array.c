/* PANDAseq -- Assemble paired FASTQ Illumina reads and strip the region between amplification primers.
     Copyright (C) 2011-2013  Andre Masella

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
#include "config.h"
#include <stdlib.h>
#include "pandaseq.h"

int panda_tweak_general_compare(
	const panda_tweak_general *a,
	const panda_tweak_general *b) {
	return a->flag - b->flag;
}

int panda_tweak_assembler_compare(
	const panda_tweak_assembler *a,
	const panda_tweak_assembler *b) {
	return a->flag - b->flag;
}

int panda_tweak_general_compare_p(
	const panda_tweak_general *const *a,
	const panda_tweak_general *const *b) {
	return panda_tweak_general_compare(*a, *b);
}

int panda_tweak_assembler_compare_p(
	const panda_tweak_assembler *const *a,
	const panda_tweak_assembler *const *b) {
	return panda_tweak_assembler_compare(*a, *b);
}

void panda_tweak_general_sort(
	const panda_tweak_general **array,
	size_t length) {
	qsort(array, length, sizeof(panda_tweak_general *), (int (*)(const void *, const void *)) panda_tweak_general_compare_p);
}

void panda_tweak_assembler_sort(
	const panda_tweak_assembler **array,
	size_t length) {
	qsort(array, length, sizeof(panda_tweak_assembler *), (int (*)(const void *, const void *)) panda_tweak_assembler_compare_p);
}

static void append_array(
	const void ***array,
	size_t *length,
	const void **additions,
	size_t additions_length) {
	size_t it;

	if (*array == NULL)
		*length = 0;

	*array = realloc(*array, (*length + additions_length) * sizeof(void *));

	for (it = 0; it < additions_length; it++) {
		(*array)[*length + it] = additions[it];
	}
	*length += additions_length;
}

void panda_tweak_general_append(
	const panda_tweak_general ***array,
	size_t *length,
	const panda_tweak_general *const *additions,
	size_t additions_length) {
	append_array((const void ***) array, length, (const void **) additions, additions_length);
	panda_tweak_general_sort(*array, *length);
}

void panda_tweak_assembler_append(
	const panda_tweak_assembler ***array,
	size_t *length,
	const panda_tweak_assembler *const *additions,
	size_t additions_length) {
	append_array((const void ***) array, length, (const void **) additions, additions_length);
	panda_tweak_assembler_sort(*array, *length);
}
