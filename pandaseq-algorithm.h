/* PANDAseq -- Assemble paired FASTQ Illumina reads and strip the region between amplification primers.
     Copyright (C) 2013  Andre Masella

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

#ifndef _PANDASEQ_ALGO_H
#        define _PANDASEQ_ALGO_H
#        ifdef __cplusplus
#                define EXTERN_C_BEGIN  extern "C" {
#                define EXTERN_C_END    }
#        else
#                define EXTERN_C_BEGIN
#                define EXTERN_C_END
#        endif
#        include <stdarg.h>
#        include <stdio.h>
#        include <stdbool.h>
EXTERN_C_BEGIN
/* === Constants === */
PANDA_EXTERN PandaAlgorithmClass *panda_algorithms;
PANDA_EXTERN size_t panda_algorithms_length;

/**
 * Add a new algorithm to the system.
 *
 * The class information must be in constant memory, as there is no way to unload an algorithm.
 */
void panda_algorithm_register(
	PandaAlgorithmClass clazz);

/* === Constructors === */
/**
 * Create a new algorithm with the specified amount of storage for the private user data.
 *
 * @clazz: the structure holding the callbacks.
 */
PandaAlgorithm panda_algorithm_new(
	PandaAlgorithmClass clazz);

/* === Getters and Setters === */

/**
 * Get the underlying algorithm class.
 */
PandaAlgorithmClass panda_algorithm_class(
	PandaAlgorithm algo);

/**
 * Get the implementation-defined data for this algorithm.
 */
void *panda_algorithm_data(
	PandaAlgorithm algo);

/* === Methods === */

/**
 * Compute the log probability of bases matching or mismatching with the
 * supplied PHRED scores.
 */
double panda_algorithm_quality_compare(
	PandaAlgorithm algorithm,
	const panda_qual *a,
	const panda_qual *b);

/**
 * Checks if an algorithm is a member of a particular class.
 */
bool panda_algorithm_is_a(
	PandaAlgorithm algo,
	PandaAlgorithmClass clazz);

#        define PANDA_ALGORITHM_IS_A(algo, clazz) panda_algorithm_is_a(algo, (PandaAlgorithmClass) ((size_t)&class_ ## clazz 1))

/**
 * Increase the reference count on an algorithm.
 */
PandaAlgorithm panda_algorithm_ref(
	PandaAlgorithm algo);

/**
 * Decrease the reference count on an algorithm.
 * @proxy: (transfer full): the algorithm to release.
 */
void panda_algorithm_unref(
	PandaAlgorithm algo);

typedef struct panda_algorithm panda_algorithm;

/* === Simple Bayes === */

PANDA_EXTERN const struct panda_algorithm_class panda_algorithm_simple_bayes_class;
#        define class_panda_algorithm_simple_bayes panda_algorithm_simple_bayes_class
typedef struct panda_algorithm panda_algorithm_simple_bayes;

/**
 * Create a simple Bayesian algorithm.
 */
PandaAlgorithm panda_algorithm_simple_bayes_new(
	void);

/**
 * The minimum error estimation in the sequence data (epsilon)
 */
double panda_algorithm_simple_bayes_get_error_estimation(
	PandaAlgorithm algorithm);
void panda_algorithm_simple_bayes_set_error_estimation(
	PandaAlgorithm algorithm,
	double q);

/* === PEAR === */

PANDA_EXTERN const struct panda_algorithm_class panda_algorithm_pear_class;
#        define class_panda_algorithm_pear panda_algorithm_pear_class
typedef struct panda_algorithm panda_algorithm_pear;

/**
 * Create a PEAR algorithm from Zhang 2013.
 */
PandaAlgorithm panda_algorithm_pear_new(
	void);

/**
 * The log probability of a random base match.
 */
double panda_algorithm_pear_get_random_base_log_p(
	PandaAlgorithm algorithm);
void panda_algorithm_pear_set_random_base_log_p(
	PandaAlgorithm algorithm,
	double log_p);

/* === RDP Maximum Likehood === */

/**
 * This algorithm uses a maximum likelihood approach to find the best overlap
 * region, and calculates an assembled Q score for each base in the overlap
 * region 
 */
PANDA_EXTERN const struct panda_algorithm_class panda_algorithm_rdp_mle_class;
#        define class_panda_algorithm_rdp_mle panda_algorithm_rdp_mle_class;
typedef struct panda_algorithm panda_algorithm_rdp_mle;

/**
 * Create a modified bayesian algorithm from RDP.
 */
PandaAlgorithm panda_algorithm_rdp_mle_new(
	void);

EXTERN_C_END
#endif
