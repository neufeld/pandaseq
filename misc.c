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
#include "config.h"
#include <unistd.h>
#if HAVE_SYS_PARAM_H
#        include <sys/param.h>
#endif
#if HAVE_SYS_SYSCTL_H
#        include <sys/sysctl.h>
#endif
#include "pandaseq.h"

char const *const panda_version(
	void) {
	return PACKAGE_STRING;
}

int panda_api_version(
	void) {
	return PANDA_API;
}

size_t panda_max_len(
	void) {
	return MAX_LEN;
}

#ifdef _WIN32
#        include <windows.h>
int panda_get_default_worker_threads(
	void) {
	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);

	return (sysinfo.dwNumberOfProcessors < 1) ? 1 : sysinfo.dwNumberOfProcessors;
}
#elif defined(_SC_NPROCESSORS_ONLN)
int panda_get_default_worker_threads(
	void) {
	int num_cpus = sysconf(_SC_NPROCESSORS_ONLN);
	return (num_cpus < 1) ? 1 : num_cpus;
}
#elif defined(HW_AVAILCPU)
#        include <unistd.h>
#        include <sys/types.h>
#        include <sys/sysctl.h>
int panda_get_default_worker_threads(
	void) {
	int num_cpus;
	int mib[4];
	size_t len = sizeof(num_cpus);
	mib[0] = CTL_HW;
	mib[1] = HW_AVAILCPU;
	sysctl(mib, 2, &num_cpu, &len, NULL, 0);

	if (num_cpu < 1) {
		mib[1] = HW_NCPU;
		sysctl(mib, 2, &numCPU, &len, NULL, 0);
	}
	return (num_cpus < 1) ? 1 : num_cpus;
}
#else
int panda_get_default_worker_threads(
	void) {
	return 1;
}
#endif
